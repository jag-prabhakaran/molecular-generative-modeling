import os.path
import json
from tensorflow import keras
import tensorflow as tf
from rdkit import Chem
from rdkit.Chem import Crippen
from rdkit.Chem import QED
from rdkit.Chem import Descriptors as descriptors
from rdkit.Chem import Lipinski
from generator import GraphWGAN, generator, discriminator, graph_to_molecule


class ModelHandler(object):
    """
    A ModelHandler defines a model with loading and inference methods that
    are used by the SageMaker serverless inference container.
    """

    def __init__(self):
        self.initialized = False
        self.model = None
        self.num_molecules = 48
        self.log_p_min = 0
        self.log_p_max = 5

    def initialize(self, context):
        """
        Initialize model. This will be called during model loading time
        :param context: Initial context contains model server system properties.
        :return:
        """
        self.initialized = True
        self.model = GraphWGAN(generator, discriminator, discriminator_steps=1)
        self.model.compile(
            optimizer_generator=keras.optimizers.Adam(5e-4),
            optimizer_discriminator=keras.optimizers.Adam(5e-4),
        )
        self.model.built = True
        model_dir = context.system_properties.get("model_dir")
        self.model.load_weights(os.path.join(model_dir, "model_weights.h5"))

    def preprocess(self, request):
        """
        Pre-process request data before feeding it to the loaded model for inference.
        :param request: JSON string of request payload.
        :return: list of preprocessed model input
        """
        print(request)
        request = json.loads(request[0]["body"])
        print(request)
        self.num_molecules = request["num_molecules"]
        self.log_p_min = request["log_p_min"]
        self.log_p_max = request["log_p_max"]
        self.qed_max = request["qed_max"]
        self.qed_min = request["qed_min"]

    def inference(self):
        """
        Internal inference method
        :return:
        """
        z = tf.random.normal((self.num_molecules, 64))
        graph = self.model.generator.predict(z)
        adjacency = tf.argmax(graph[0], axis=1)
        adjacency = tf.one_hot(adjacency, depth=4, axis=1)
        adjacency = tf.linalg.set_diag(adjacency, tf.zeros(tf.shape(adjacency)[:-1]))
        features = tf.argmax(graph[1], axis=2)
        features = tf.one_hot(features, depth=5, axis=2)
        return [
            graph_to_molecule([adjacency[i].numpy(), features[i].numpy()])
            for i in range(self.num_molecules)
        ]

    def postprocess(self, model_output):
        """
        Post-processing function to apply to the model output before returning to client
        :param model_output: Output from the inference method
        :return: Output after post-processing step
        """
        smiles_list = [
            {
                "smile": Chem.MolToSmiles(mol),
                "logP": Crippen.MolLogP(mol),
                "qed": QED.default(mol),
                "mol_weight": descriptors.ExactMolWt(mol),
                "num_h_donors": Lipinski.NumHDonors(mol),
            }
            for mol in model_output
            if mol is not None
        ]
        filtered_smiles = list(
            filter(
                lambda x: x["logP"] >= self.log_p_min
                and x["logP"] <= self.log_p_max
                and x["qed"] >= self.qed_min
                and x["qed"] <= self.qed_max,
                smiles_list,
            )
        )
        # return as list to keep sagemaker mms happy
        return [json.dumps({"smiles": smiles_list, "filtered_smiles": filtered_smiles})]

    def ping(self):
        """
        Ping to get system health
        :return:
        """
        # TODO: Implement actual health check
        return "PONG"

    def handle(self, data, context):
        """
        Call preprocess, inference and post-process functions
        :param data:
        :param context:
        :return:
        """
        self.preprocess(data)
        raw_output = self.inference()
        return self.postprocess(raw_output)


_service = ModelHandler()


def handle(data, context):
    """
    Sagemaker inference handler function
    :param data:
    :param context:
    :return:
    """
    if not _service.initialized:
        _service.initialize(context)

    if data is None:
        return None

    return _service.handle(data, context)
