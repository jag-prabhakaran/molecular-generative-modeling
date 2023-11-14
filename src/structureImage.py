import json
from rdkit import Chem
from rdkit.Chem import Draw

def generate_image(smiles, output_path):
    mol = Chem.MolFromSmiles(smiles)
    img = Draw.MolToImage(mol, size=(300, 300))
    img.save(output_path)

def process_json(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)

    for entry in data:
        smiles = entry.get('smiles')
        image_id = entry.get('id')
        output_path = f"path/to/generated/images/{image_id}.png"
        generate_image(smiles, output_path)

if __name__ == "__main__":
    json_file_path = "path/to/your/sample-response.json"
    process_json()
