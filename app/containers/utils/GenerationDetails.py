from typing import Union
from collections import defaultdict
import json
from math import inf
from rdkit import Chem
from rdkit.Chem import Mol, Crippen, QED, Descriptors, Lipinski

PROP_TO_FUNCTION_MAP = {
    "logP": Crippen.MolLogP,
    "qed": QED.qed,
    "mol_weight": Descriptors.ExactMolWt,
    "num_h_donors": Lipinski.NumHDonors,
}


class GenerationDetails(object):
    """
    Class representing the details of molecule generation.

    Attributes:
        prop_filter_ranges (dict[str, tuple[float | str]]): A dictionary mapping property names to their filter ranges.
        input_smile (str): The input SMILE string.
        output_molecules (dict[str, dict[str, float | str]]): A dictionary mapping SMILE strings to property dictionaries.
        canonical_output_smiles (list[str]): A list of canonical SMILE strings.
        generation_upper_bound (int): The upper bound on the number of generated molecules.

    Methods:
        __init__(self, input_json: str, filter_props: list[str]) -> None: Initializes the GenerationDetails object.
        filter_output(self, raw_model_smiles: list[str]): Filters the raw model output based on property ranges.
        generate_json(self) -> str: Generates a JSON string representation of the output molecules.
    """

    prop_filter_ranges: dict[str, tuple[Union[float, str]]]
    input_smile: str
    output_molecules: dict[str, dict[str, Union[float, str]]]
    canonical_output_smiles: list[str]
    generation_upper_bound: int

    def __init__(self, input_json: str, filter_props: list[str]) -> None:
        """
        Initializes the GenerationDetails object.

        Args:
            input_json (str): The input JSON string containing the input SMILE and filter properties.
            filter_props (list[str]): A list of property names to filter the output molecules.

        Returns:
            None
        """
        self.input_smile = input_json["input_smile"]
        self.generation_upper_bound = input_json["num_molecules"]
        self.prop_filter_ranges = {}
        self.output_molecules = defaultdict(dict)
        for prop in filter_props:
            if prop in input_json.keys():
                self.prop_filter_ranges[prop] = input_json[prop]
            else:
                self.prop_filter_ranges[prop] = [-inf, inf]

    def filter_output(self, raw_model_smiles: list[str]):
        """
        Filters the output molecules based on property ranges.

        Args:
            raw_model_smiles (list[str]): A list of raw model-generated SMILE strings.

        Returns:
            None
        """
        for smile in raw_model_smiles:
            print(smile)
            if smile not in self.output_molecules.keys():
                mol = Chem.MolFromSmiles(smile)
                # check if SMILE is valid
                if mol:
                    prop_dict = {}
                    valid = True
                    for prop, prop_range in self.prop_filter_ranges.items():
                        value = PROP_TO_FUNCTION_MAP[prop](mol)
                        prop_dict[prop] = value
                        if not prop_range[0] <= value <= prop_range[1]:
                            valid = False
                    if valid:
                        self.output_molecules[smile] = prop_dict

    def generate_json(self) -> str:
        """
        Generates a JSON string representation of the output molecules.

        Returns:
            str: The JSON string representation of the output molecules.
        """
        return json.dumps(self.output_molecules)
