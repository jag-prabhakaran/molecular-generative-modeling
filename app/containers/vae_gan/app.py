from flask import Flask, jsonify, request
from rdkit.Chem import Crippen

import tensorflow as tf
from generator import *

app = Flask(__name__)

new_model = GraphWGAN(generator, discriminator, discriminator_steps=1)
new_model.compile(
    optimizer_generator=keras.optimizers.Adam(5e-4),
    optimizer_discriminator=keras.optimizers.Adam(5e-4),
)
new_model.built = True
savedModel = new_model.load_weights('model_weights.h5')

def sample(generator, batch_size):
    z = tf.random.normal((batch_size, LATENT_DIM))
    graph = generator.predict(z)
    adjacency = tf.argmax(graph[0], axis=1)
    adjacency = tf.one_hot(adjacency, depth=BOND_DIM, axis=1)
    adjacency = tf.linalg.set_diag(adjacency, tf.zeros(tf.shape(adjacency)[:-1]))
    features = tf.argmax(graph[1], axis=2)
    features = tf.one_hot(features, depth=ATOM_DIM, axis=2)
    return [
        graph_to_molecule([adjacency[i].numpy(), features[i].numpy()])
        for i in range(batch_size)
    ]

@app.route('/generate_molecules', methods=['POST'])
def generate_molecules():
    data = request.get_json()
    molecules = sample(new_model.generator, batch_size=48)
    smiles_list = [[Chem.MolToSmiles(mol), Crippen.MolLogP(mol)] for mol in molecules if mol is not None]
    filtered_smiles = [x for x in smiles_list if data["log_p_min"] <= x[1] <= data["log_p_max"]]
    return jsonify({
        'smiles': smiles_list,
        'filtered_smiles': filtered_smiles
    })

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0')
