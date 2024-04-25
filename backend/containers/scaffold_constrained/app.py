from flask import Flask, request, jsonify
from rdkit import Chem
from rdkit.Chem import Crippen
from model_files.scaffold_constrained_model import scaffold_constrained_RNN
from model_files.data_structs import Vocabulary
from model_files.utils import seq_to_smiles
import torch

voc = Vocabulary(init_from_file="data/DistributionLearningBenchmark/Voc")
Agent = scaffold_constrained_RNN(voc)
Agent.rnn.load_state_dict(
    torch.load(
        "data/DistributionLearningBenchmark/Prior_ChEMBL_randomized.ckpt",
        map_location=lambda storage, loc: storage,
    )
)
# Create a Flask web application
app = Flask(__name__)

print("Model loaded")


@app.route("/run_model_inference", methods=["POST"])
def run_model_inference():
    data = request.get_json()
    seqs, agent_likelihood, entropy = Agent.sample(
        pattern=data["scaffold_smile"], batch_size=50
    )
    smiles = seq_to_smiles(seqs, voc)
    mols = []
    for smile in smiles:
        mol = Chem.MolFromSmiles(smile)
        if mol:
            mols.append(mol)
    mol_smiles = [(Chem.MolToSmiles(mol), Crippen.MolLogP(mol)) for mol in mols]
    filtered_mol_smiles = [
        x for x in mol_smiles if data["log_p_min"] <= x[1] <= data["log_p_max"]
    ]
    return jsonify({"smiles": mol_smiles, "filtered_smiles": filtered_mol_smiles}), 200


# Run the application on a specific host and port
if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000)
