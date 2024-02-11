interface Molecule {
    get_substruct_matches(qmol: Molecule): string
    "smile": string,
    "logP": number,
    "qed": number,
    "mol_weight": number,
    "num_h_donors": number
}  

export default Molecule