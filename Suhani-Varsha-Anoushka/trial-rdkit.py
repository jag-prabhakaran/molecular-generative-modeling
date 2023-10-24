from rdkit import Chem
from rdkit.Chem import Draw

smiles_string = "N[C@@H](CCC(=O)N[C@@H](CS)C(=O)NCC(=O)O)C(=O)O"
molecule = Chem.MolFromSmiles(smiles_string)
img = Draw.MolToImage(molecule)

img.save("d:\\University Studies\\Fall 2023\\Data Mine - Merck\\molecular-generative-modeling\\Suhani-Varsha-Anoushka\\molecule.png")