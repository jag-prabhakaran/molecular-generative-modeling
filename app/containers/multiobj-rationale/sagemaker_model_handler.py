import os.path
import json
from rdkit import Chem
from rdkit.Chem import Crippen
from rdkit.Chem import QED
from rdkit.Chem import Descriptors as descriptors
from rdkit.Chem import Lipinski
import torch
import random
import sys

from torch.utils.data import DataLoader

from fuseprop import *

from types import SimpleNamespace


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
        self.args = None
        self.qed_max = 1
        self.qed_min = 0
        self.rationale = []

    def initialize(self, context):
        """
        Initialize model. This will be called during model loading time
        :param context: Initial context contains model server system properties.
        :return:
        """
        self.initialized = True
        args = SimpleNamespace(
            atom_vocab=common_atom_vocab,
            model="ckpt/chembl-h400beta0.3/model.20",
            num_decode=100,
            seed=1,
            rnn_type="LSTM",
            hidden_size=400,
            embed_size=400,
            batch_size=20,
            latent_size=20,
            depth=10,
            diter=3,
        )
        self.args = args
        random.seed(1)
        self.model = AtomVGNN(args).cpu()
        model_dir = context.system_properties.get("model_dir")
        model_ckpt = torch.load(
            os.path.join(model_dir, args.model), map_location=torch.device("cpu")
        )
        print("loading pre-trained model", file=sys.stderr)
        self.model.load_state_dict(model_ckpt)
        self.model.eval()
        torch.manual_seed(args.seed)

    def preprocess(self, request):
        """
        Pre-process request data before feeding it to the loaded model for inference.
        :param request: JSON string of request payload.
        :return: list of preprocessed model input
        """
        request = json.loads(request[0]["body"])
        self.rationale = request["rationale"]
        self.num_molecules = request["num_molecules"]
        self.args.num_decode = self.num_molecules
        self.log_p_min = request["log_p_min"]
        self.log_p_max = request["log_p_max"]
        self.qed_max = request["qed_max"]
        self.qed_min = request["qed_min"]

    def inference(self):
        """
        Internal inference method
        :return:
        """
        dataset = SubgraphDataset(
            self.rationale,
            self.args.atom_vocab,
            self.args.batch_size,
            self.args.num_decode,
        )
        loader = DataLoader(
            dataset,
            batch_size=1,
            shuffle=False,
            num_workers=0,
            collate_fn=lambda x: x[0],
        )
        with torch.no_grad():
            return [self.model.decode(init_smiles) for init_smiles in loader]

    def postprocess(self, model_output):
        """
        Post-processing function to apply to the model output before returning to client
        :param model_output: Output from the inference method
        :return: Output after post-processing step
        """
        input_smiles = self.rationale
        output_objects = [
            {
                "input_smile": input_smile,
                "output_smile": output_smile,
                "logP": Crippen.MolLogP(Chem.MolFromSmiles(output_smile)),
                "qed": QED.default(Chem.MolFromSmiles(output_smile)),
                "mol_weight": descriptors.ExactMolWt(Chem.MolFromSmiles(output_smile)),
                "num_h_donors": Lipinski.NumHDonors(Chem.MolFromSmiles(output_smile)),
            }
            for input_smile, output_smiles in zip(input_smiles, model_output)
            for output_smile in output_smiles
        ]
        filtered_output_objects = list(
            filter(
                lambda x: x["log_p"] >= self.log_p_min
                and x["log_p"] <= self.log_p_max
                and x["qed"] >= self.qed_min
                and x["qed"] <= self.qed_max,
                output_objects,
            )
        )

        # return as list to keep sagemaker mms happy
        return [
            json.dumps(
                {
                    "output_objects": output_objects,
                    "filtered_output_objects": filtered_output_objects,
                }
            )
        ]

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
