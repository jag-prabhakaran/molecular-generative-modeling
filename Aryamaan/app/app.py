import hashlib
from flask import Flask, request, jsonify
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy.dialects.postgresql import TEXT, REAL, UUID, ARRAY
from sqlalchemy.orm import Mapped, mapped_column
from sqlalchemy import select
from rdkit import Chem
from rdkit.Chem import Crippen
from rdkit.Chem import AllChem, Draw
from rdkit import DataStructs
from scaffold_constrained_model import scaffold_constrained_RNN
from data_structs import Vocabulary, Experience
from utils import Variable, seq_to_smiles, fraction_valid_smiles, unique
import torch
import os
from dotenv import load_dotenv

load_dotenv()

# Create a Flask web application
app = Flask(__name__)

# Configure the database connection (replace with your own PostgreSQL URI)
app.config[
    'SQLALCHEMY_DATABASE_URI'] = os.getenv('POSTGRES_URL_SECRET')
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False  # Disable tracking modifications for SQLAlchemy

# Create a SQLAlchemy database object
db = SQLAlchemy(app)

voc = Vocabulary(init_from_file="data/DistributionLearningBenchmark/Voc")
Agent = scaffold_constrained_RNN(voc)
Agent.rnn.load_state_dict(torch.load('data/DistributionLearningBenchmark/Prior_ChEMBL_randomized.ckpt', map_location=lambda storage, loc: storage))
print("Model loaded")


# Define models for the tables in the database
class GeneratedMolecule(db.Model):
    smile: Mapped[str] = mapped_column(TEXT, primary_key=True)
    predicted_log_p: Mapped[float] = mapped_column(REAL, nullable=False)

    def __init__(self, smile: Mapped[str], predicted_log_p: Mapped[float]):
        self.predicted_log_p = predicted_log_p
        self.smile = smile

    def as_dict(self):
        return {c.name: getattr(self, c.name) for c in self.__table__.columns}


class ModelInference(db.Model):
    id: Mapped[str] = mapped_column(UUID, primary_key=True)
    generated_smile_array: Mapped[list] = mapped_column(ARRAY(TEXT), nullable=False)

    def __init__(self, id: Mapped[str], generated_smile_array: Mapped[list]):
        self.id = id
        self.generated_smile_array = generated_smile_array

    def as_dict(self):
        return {c.name: getattr(self, c.name) for c in self.__table__.columns}


class ModelInput(db.Model):
    id: Mapped[str] = mapped_column(UUID, primary_key=True)
    scaffold_smile: Mapped[str] = mapped_column(TEXT, nullable=False)
    log_p_min: Mapped[float] = mapped_column(REAL, nullable=False)
    log_p_max: Mapped[float] = mapped_column(REAL, nullable=False)

    def __init__(self, id, scaffold_smile, log_p_min, log_p_max):
        self.log_p_max = log_p_max
        self.log_p_min = log_p_min
        self.scaffold_smile = scaffold_smile
        self.id = id

    def as_dict(self):
        return {c.name: getattr(self, c.name) for c in self.__table__.columns}


# Define a route to insert data into the database using a POST request
@app.route('/insert_model_input', methods=['POST'])
def insert_data():
    data = request.get_json()
    if "scaffold_smile" not in data or "log_p_min" not in data or "log_p_max" not in data:
        return jsonify("Bad Request! Malformed data."), 400
    else:
        scaffold_smile = data['scaffold_smile']
        log_p_min = data['log_p_min']
        log_p_max = data['log_p_max']
        id = hashlib.md5((scaffold_smile + str(log_p_min) + str(log_p_max)).encode('utf-8')).hexdigest()
        input_to_insert = ModelInput(id, scaffold_smile, log_p_min, log_p_max)
        db.session.add(input_to_insert)
        db.session.commit()
        return jsonify(("Data inserted successfully!", {
            "id": id,
        })), 200

@app.route('/list_model_inputs', methods=['GET'])
def list_model_inputs():
    results = db.session.scalars(select(ModelInput))
    return [x.as_dict() for x in results.all()]

@app.route('/list_model_inferences', methods=['GET'])
def list_model_inferences():
    results = db.session.scalars(select(ModelInference))
    return [x.as_dict() for x in results.all()]

@app.route('/list_generated_molecules', methods=['GET'])
def list_generated_molecules():
    results = db.session.scalars(select(GeneratedMolecule))
    return [x.as_dict() for x in results.all()]

@app.route('/insert_model_inference', methods=['POST'])
def insert_model_inference():
    data = request.get_json()
    if "generated_smile_array" not in data or "id" not in data:
        return jsonify("Bad Request! Malformed data."), 400
    else:
        generated_smile_array = data['generated_smile_array']
        id = data['id']
        inference_to_insert = ModelInference(id, generated_smile_array)
        db.session.add(inference_to_insert)
        db.session.commit()
        return jsonify("Data inserted successfully"), 200

@app.route('/insert_generated_molecule', methods=['POST'])
def insert_generated_molecule():
    data = request.get_json()
    if "smile" not in data or "predicted_log_p" not in data:
        return jsonify("Bad Request! Malformed data."), 400
    else:
        smile = data['smile']
        predicted_log_p = data['predicted_log_p']
        molecule_to_insert = GeneratedMolecule(smile, predicted_log_p)
        db.session.add(molecule_to_insert)
        db.session.commit()
        return jsonify("Data inserted successfully"), 200

@app.route('/run_model_inference', methods=['POST'])
def run_model_inference():
    data = request.get_json()
    if "scaffold_smile" not in data or "log_p_min" not in data or "log_p_max" not in data:
        return jsonify("Bad Request! Malformed data."), 400
    else:
        seqs, agent_likelihood, entropy = Agent.sample(pattern=data["scaffold_smile"], batch_size=50)
        smiles = seq_to_smiles(seqs, voc)
        mols = []
        for smile in smiles:
            mol = Chem.MolFromSmiles(smile)
            if mol:
                mols.append(mol)
        mol_smiles = [(Chem.MolToSmiles(mol), Crippen.MolLogP(mol)) for mol in mols]
        filtered_mol_smiles = [x for x in mol_smiles if data["log_p_min"] <= x[1] <= data["log_p_max"]]
        return jsonify({
            "mol_smiles": mol_smiles,
            "filtered_mol_smiles": filtered_mol_smiles
        }), 200


# Run the application on a specific host and port
if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)
    with app.app_context():
        db.create_all()
