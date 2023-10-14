from flask import Flask, jsonify

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

@app.route('/generate_molecules', methods=['GET'])
def generate_molecules():
    molecules = sample(new_model.generator, batch_size=48)
    smiles_list = [Chem.MolToSmiles(mol) for mol in molecules if mol is not None]
    return jsonify({'generated_smiles': smiles_list})

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0')
