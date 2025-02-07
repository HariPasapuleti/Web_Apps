# Molecular Solubility Prediction Web App

This repository contains the code for a **Molecular Solubility Prediction Web App**, which predicts the **Solubility (LogS)** values of molecules based on their SMILES representation. The app utilizes machine learning and molecular descriptors to make predictions and visualize the results.

## Features

- **Molecular Descriptor Calculation**: The app computes key molecular descriptors such as MolLogP, MolWt, and AromaticProportion.
- **SMILES Input**: Users can input molecular structures as SMILES strings.
- **Solubility Prediction**: The app uses a pre-trained model to predict the solubility (LogS) of the molecules.
- **Interactive Interface**: Built with Streamlit for easy use and interaction.

## Technologies Used

- **Streamlit**: For building the web application.
- **RDKit**: For calculating molecular descriptors from SMILES strings.
- **scikit-learn**: For loading the pre-trained machine learning model to make predictions.
- **Pandas & NumPy**: For data manipulation.
- **Pickle**: For saving and loading the trained model.

## Getting Started

To get started with this project locally, follow the steps below:

### Prerequisites

Make sure you have the following installed:

- Python 3.7+
- pip (Python package manager)

### Installation

1. Clone the repository to your local machine:

    ```bash
    git clone https://github.com/HariPasapuleti/Web_Apps.git
    cd Web_Apps
    ```

2. Install the necessary dependencies:

    ```bash
    pip install -r requirements.txt
    ```

3. Run the app:

    ```bash
    streamlit run solubility_prediction_app.py
    ```

4. Open the app in your browser.

### Input

- Users can input molecular structures in the form of SMILES strings through the sidebar. The default SMILES input is:
NCCCC CCC CN


### Output

The app will display the following outputs:

- **Molecular Descriptors**: LogP, Molecular Weight, Number of Rotatable Bonds, and Aromatic Proportion.
- **Predicted LogS Value**: The predicted solubility value for the input molecule.

## Model Details

The app uses a pre-trained machine learning model saved in the file `solubility_model.pkl`. The model takes in molecular descriptors and predicts the solubility of the molecule (LogS).

## Deployed Link

You can access the deployed web app here:
[Deployed Molecular Solubility Prediction Web App](https://molecular-solubility.streamlit.app/)
