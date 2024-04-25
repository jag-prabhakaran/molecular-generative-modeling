import os.path
import json
from rdkit import Chem
from rdkit.Chem import Crippen
from rdkit.Chem import QED
from rdkit.Chem import Descriptors as descriptors
from rdkit.Chem import Lipinski
from model_files.scaffold_constrained_model import scaffold_constrained_RNN
from model_files.data_structs import Vocabulary
import torch
from model_files.utils import seq_to_smiles
from utils import GenerationDetails


class ModelHandler(object):
    """
    A ModelHandler defines a model with loading and inference methods that
    are used by the SageMaker serverless inference container.
    """

    def __init__(self):
        self.initialized = False
        self.model = None
        self.gen_details = None

    def initialize(self, context):
        """
        Initialize model. This will be called during model loading time
        :param context: Initial context contains model server system properties.
        :return:
        """
        self.initialized = True
        model_dir = context.system_properties.get("model_dir")
        voc = Vocabulary(
            init_from_file=os.path.join(
                model_dir, "data/DistributionLearningBenchmark/Voc"
            )
        )
        self.model = scaffold_constrained_RNN(voc)
        self.model.rnn.load_state_dict(
            torch.load(
                os.path.join(
                    model_dir,
                    "data/DistributionLearningBenchmark/Prior_ChEMBL_randomized.ckpt",
                ),
                map_location=lambda storage, loc: storage,
            )
        )

    def preprocess(self, request):
        """
        Pre-process request data before feeding it to the loaded model for inference.
        :param request: JSON string of request payload.
        :return: list of preprocessed model input
        """
        request = json.loads(request[0]["body"])
        self.gen_details = GenerationDetails(
            input_json=request,
            filter_props=["logP", "qed", "mol_weight", "num_h_donors"],
        )

    def inference(self):
        """
        Internal inference method
        :return:
        """
        seqs, agent_likelihood, entropy = self.model.sample(
            pattern=self.gen_details.input_smile,
            batch_size=self.gen_details.generation_upper_bound,
        )
        return seqs

    def postprocess(self, model_output):
        """
        Post-processing function to apply to the model output before returning to client
        :param model_output: Output from the inference method
        :return: Output after post-processing step
        """
        smiles = seq_to_smiles(model_output, self.model.voc)
        self.gen_details.filter_output(smiles)
        # return as list to keep sagemaker mms happy
        return [self.gen_details.generate_json()]

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
