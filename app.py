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
  aa_count = []
  for i in aromatic_atoms:
    if i==True:
      aa_count.append(1)
  AromaticAtom = sum(aa_count)
  HeavyAtom = Descriptors.HeavyAtomCount(m)
  AR = AromaticAtom/HeavyAtom
  return AR

def generate(smiles, verbose=False): # Verbose provides the extra information about its internal process {True or False}

    moldata= []
    for elem in smiles:
        mol=Chem.MolFromSmiles(elem) #Convert the elem in smiles string into RDKit mol object
        moldata.append(mol)

    baseData= np.arange(1,1)
    #Creates an empty array because the start and stop values are the same, resulting in no values to include in the array.

    i=0
    for mol in moldata:

        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_MolWt = Descriptors.MolWt(mol)
        desc_NumRotatableBonds = Descriptors.NumRotatableBonds(mol)
        desc_AromaticProportion = AromaticProportion(mol)

        row = np.array([desc_MolLogP,
                        desc_MolWt,
                        desc_NumRotatableBonds,
                        desc_AromaticProportion])

        if(i==0):
            baseData=row
        else:
            baseData=np.vstack([baseData, row])
            # Uses np.vstack to vertically stack baseData with the new row. 
            # np.vstack takes a list of arrays and stacks them on top of each other, creating a new array where each row is added as a new row in baseData.
        i=i+1

    columnNames=["MolLogP","MolWt","NumRotatableBonds","AromaticProportion"]
    descriptors = pd.DataFrame(data=baseData,columns=columnNames)

    return descriptors


image = Image.open('Solubility_Prediction.png')

st.image(image, use_container_width=True)

st.write("""
# Molecular Solubility Prediction Web App

This app predicts the **Solubility (LogS)** values of molecules!

***
""")

# Input molecules (Side Panel)

st.sidebar.header('User Input Features')

SMILES_input = "NCCCC\nCCC\nCN" #Butylamine\Propane\Methylamaine

SMILES = st.sidebar.text_area("SMILES input", SMILES_input)
SMILES = "C\n" + SMILES #Adds C as a dummy, first item
SMILES = SMILES.split('\n')

st.header('Input SMILES')
SMILES[1:] # Skips the dummy first item


st.header('Computed molecular descriptors')
X = generate(SMILES)
X[1:] # Skips the dummy first item

load_model = pickle.load(open('solubility_model.pkl', 'rb')) # rb Binary read mode
 
prediction = load_model.predict(X)
#prediction_proba = load_model.predict_proba(X)

st.header('Predicted LogS values')
prediction[1:] # Skips the dummy first item
