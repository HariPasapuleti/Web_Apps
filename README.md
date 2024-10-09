# Execution
For execution run `streamlit run app.py`

## Aromatic proportion

Aromatic proportion refers to the fraction of aromatic atoms in a molecule relative to the total number of atoms. 
* Aromatic Proportion = Total Number of Atoms / Number of Aromatic Atoms.
​

# Points
1. `pd.DataFrame():` This is a constructor function from the Pandas library that creates a DataFrame, which is a two-dimensional, size-mutable, and potentially heterogeneous tabular data structure.

2. `data=baseData:` baseData: This is the data that will populate the DataFrame. It should be a 2D array-like structure (e.g., a list of lists, a NumPy array, or another DataFrame). Each sub-array or row in baseData represents a row in the DataFrame.
columns=columnNames:

3. `columnNames:` This is a list of column names that will be assigned to the DataFrame’s columns. It should have the same length as the number of columns in baseData.


4. `st.image():` Function to display images in a Streamlit app.
5. `image:` The image to be displayed, which can be a file path, URL, or array.
6. `use_column_width=True:` 
* Ensures the image width adjusts to the width of the column in the Streamlit app.
* This functionality helps in creating visually appealing and responsive applications by automatically adjusting the size of images to fit the available space.


7. `st.sidebar.text_area("SMILES input", SMILES_input):` 
* Function: Creates a text area in the sidebar for user input.
* Purpose: Allows users to input or edit SMILES strings. The default value is set to the SMILES_input variable.
8. `SMILES = "C\n" + SMILES:` 
* Action: Prepends "C\n" to the existing SMILES input.
* Purpose: Adds a dummy SMILES string, "C", to the beginning of the input. This is done to ensure that the resulting list of SMILES strings starts with a valid entry.
* Adding "C": Ensures that the list of SMILES strings starts with a valid entry, which can prevent errors and handle edge cases in data processing.
9. `SMILES = SMILES.split('\n'):` 
* Function: Splits the SMILES string into a list of strings using the newline character (\n) as the delimiter.
* Purpose: Converts the multi-line input into a list of SMILES strings.


10. Open the File: `open('solubility_model.pkl', 'rb')` opens the specified file in binary read mode. The file should be in the same directory as your script or notebook unless a path is specified.
11. Load the Model: `pickle.load()` reads the content of the file and deserializes it. This requires that the file solubility_model.pkl was previously saved using pickle.dump().
12. Assign to Variable: The deserialized model is assigned to the variable load_model, which now holds the model object that was trained and saved previously.

Example1
```python
import pandas as pd
import numpy as np

# Example data
baseData = np.array([
    [1, 2, 3],
    [4, 5, 6],
    [7, 8, 9]
])

# Column names
columnNames = ['Column1', 'Column2', 'Column3']

# Create the DataFrame
descriptors = pd.DataFrame(data=baseData, columns=columnNames)

print(descriptors)
```
output:
```
   Column1  Column2  Column3
0        1        2        3
1        4        5        6
2        7        8        9
```


Example2: 
If a user inputs the following in the text area:

```
NCCCC
CCC
CN
```
With the dummy "C" added, the processed result will be:

`SMILES = ["C", "NCCCC", "CCC", "CN"]`
This ensures that the list starts with "C", followed by the user-provided SMILES strings.

