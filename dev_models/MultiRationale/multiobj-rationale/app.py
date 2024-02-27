import streamlit as st
from streamlit_ketcher import st_ketcher
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw
from decode import *
from rdkit.Chem.Descriptors import qed
import re

title = st.header("Draw the Molecule")
input_smiles = st_ketcher()
# st.text(decode(smiles,3))
# smiles = "OCc1cc[c:1]c(-c2ncccn2)c1"
p = re.compile(r"([a-zA-z][0-9]*)\(\*\)")
smiles = p.sub("[\\1:1]", input_smiles)
st.write(input_smiles, smiles)


n_mol = 3

with st.sidebar:
    n_mol = st.number_input(
        "Number of Molecules", value=5, placeholder="Type a number..."
    )
    genre = st.radio(
        "Property",
        ["QED", "SA", "LogP"],
        captions=[
            "Quantitative Estimate of Druglikeness",
            "Synthetic Accessibility",
            "Partition Coefficient",
        ],
    )


def calc_qed(smiles_list):
    scores = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            scores.append(0)
        else:
            scores.append(qed(mol))
    return np.float32(scores)


def v_spacer(height, sb=False) -> None:
    for _ in range(height):
        if sb:
            st.sidebar.write("\n")
        else:
            st.write("\n")


col1, col2 = st.columns(2)
with col1:
    st.header("Input Molecule")

with col2:
    st.header("QED Score")

with col1:
    m = Chem.MolFromSmiles(input_smiles)
    im = Draw.MolToImage(m, size=(200, 200), fitImage=True, useSVG=True)
    st.image(im)
with col2:
    v_spacer(height=5)
    st.markdown("""### {}""".format(qed(m)))
    v_spacer(height=5)


if st.button("Generate"):
    # title = st.header('Output Molecules')
    output_mols = decode(smiles, n_mol)
    qed_scores = calc_qed(output_mols)
    indices = np.argsort(qed_scores)[::-1]

    with col1:
        st.header("Output Molecule")

    with col2:
        st.header("QED Score")

    for indx in indices:
        with col1:
            m = Chem.MolFromSmiles(output_mols[indx])
            im = Draw.MolToImage(m, size=(200, 200), fitImage=True, useSVG=True)
            st.image(im)
        with col2:
            v_spacer(height=5)
            st.markdown("""### {}""".format(qed_scores[indx]))
            v_spacer(height=5)
