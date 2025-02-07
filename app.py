import numpy as np
import pandas as pd
import streamlit as st
import pickle
from PIL import Image
from rdkit import Chem
from rdkit.Chem import Descriptors


## Calculate molecular descriptors
def AromaticProportion(m):
    aromatic_atoms = [m.GetAtomWithIdx(i).GetIsAromatic() for i in range(m.GetNumAtoms())]
    AromaticAtom = sum(1 for i in aromatic_atoms if i)
    HeavyAtom = Descriptors.HeavyAtomCount(m)
    AR = AromaticAtom / HeavyAtom if HeavyAtom > 0 else 0
    return AR

def generate(smiles):
    moldata = [Chem.MolFromSmiles(elem) for elem in smiles]
    baseData = []

    for mol in moldata:
        if mol:  # Check for valid molecules
            row = [
                Descriptors.MolLogP(mol),
                Descriptors.MolWt(mol),
                Descriptors.NumRotatableBonds(mol),
                AromaticProportion(mol)
            ]
            baseData.append(row)
    
    columnNames = ["MolLogP", "MolWt", "NumRotatableBonds", "AromaticProportion"]
    return pd.DataFrame(data=baseData, columns=columnNames)


# Load the app banner image
image = Image.open('Solubility_Prediction.png')
st.image(image)

st.write("""
# Molecular Solubility Prediction Web App

This app predicts the **Solubility (LogS)** values of molecules!

***
""")

# Sidebar input
st.sidebar.header('User Input Features')
SMILES_input = "NCCCC\nCCC\nCN"  # Default example
SMILES = st.sidebar.text_area("SMILES input", SMILES_input).split('\n')

st.header('Input SMILES')
st.write(SMILES)

# Compute molecular descriptors
st.header('Computed Molecular Descriptors')
X = generate(SMILES)
st.write(X)

# Load the trained model
try:
    load_model = pickle.load(open('solubility_model.pkl', 'rb'))
    # Dummy prediction for validation
    load_model.predict([[0, 0, 0, 0]])
except Exception as e:
    st.error(f"Model validation failed: {e}")
    st.stop()

# Make predictions
try:
    predictions = load_model.predict(X.to_numpy())
    st.header('Predicted LogS values')
    st.write(predictions)
except Exception as e:
    st.error(f"Prediction failed: {e}")
