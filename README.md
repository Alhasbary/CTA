# Compound-Target Activity (CTA) Prediction Tool

The Compound-Target Activity (CTA) Prediction Tool is an open-source ligand-based program designed to identify target information associated with chemical compounds. This command-line tool, written in Python and built upon RDKit, provides a comprehensive suite of tools for preparing and managing compound databases.

## Key Features

### 1. **extractDB: Extraction of High-Quality Data from ChEMBL**

This tool allows you to retrieve high-quality bioactivity data from the ChEMBL database. It comes with predefined selection constraints and offers the flexibility to apply user-specified constraints, such as confidence scores and activity values for assays.

### 2. **processMS: Processing of Molecular Structures**

Responsible for cleaning SMILES strings using RDKit and MolVS, this tool performs essential preprocessing steps, including structure normalization, desalting, and charge neutralization. Users have the option to disable desalting and charge neutralization if needed.

### 3. **conductSS: Conducting Chemical Similarity Searching**

This tool enables chemical similarity searches, a fundamental step in ligand-based CTA prediction. It retrieves potential targets for each given SMILES string. Users can customize fingerprint settings, including type, nBits, and radius parameters that control substructure size. A similarity threshold (Tc) can also be specified.

### 4. **createCTA: Creating Compound-Target Activity Pair Dataset**

Creating a compound-target activity pair dataset is crucial for developing ligand-based CTA prediction models, including machine learning models capable of identifying target information for chemical compounds.

## System Requirements

- Python 3
- General Python Modules: numpy, pandas, joblib
- Special Python Modules: RDKit and MolVS
- ChEMBL SQLite database

## Installation

1. Download the repository and extract files to your desired location.
2. Download the ChEMBL SQLite dataset from [ChEMBL Downloads](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/).
3. Open a command-line window (Windows) or a terminal (Linux).
4. Run the Python scripts as described below.

## Usage

### Create a Custom mini-ChEMBL SQLite Database

Create a custom mini-ChEMBL SQLite database tailored to fulfill specific application requirements. In our case, it significantly reduces the storage size from 22.4 GB (ChEMBL32) to just 714 MB.

- Inputs:
   - data: Full path to the ChEMBL dataset [Required].
   - destination: Full path to save the custom ChEMBL version as 'mini_chembl.db' [Optional].

- Output:
   - 'mini_chembl.db' is stored in the specified output directory, containing customized information for the application's use.

- Usage:
   - For help:
     ```
     python mini_chembl.py -h
     ```
   - To run the program:
     ```
     python mini_chembl.py --data=FullpathToDatabaseFile --destination=FullPathToSavedDatabase
     ```

- Usage Example:
    1. Extract the dataset: `tar -xvzf chembl_32_sqlite.tar.gz`
    2. Run the program:
     ```
     python mini_chembl.py --data=chembl_32/chembl_32_sqlite/chembl_32.db --destination=input
     ```

### Run Compound-Target Activity (CTA) Prediction Program

- Inputs:
   - List of optional parameters.
   - input: Full path to the data folder containing SMILES string lists (e.g., `*_smiles.csv`) and the mini-ChEMBL database [Optional]. It may contain more than one list. The tool processes them sequentially. Each file must be a CSV file that contains smiles strings in a column named 'smiles' and compound id in a column named 'smiles_id'.
   - output: Full path to save the results [Optional].

- Output:
   - `*_potential_targets_ss.csv`: Contains all likely targets for the chemical compounds in each dataset.
   - `*_CTA_dataset.csv`: Stores the CTA dataset.

- Usage:
   - For help:
     ```
     python main.py -h
     ```
   - To run the program:
     ```
     python main.py [parameters] --input=FullpathToInputFiles --output=FullPathToSaveResults
     ```

- Usage Example:
    ```
    python main.py
    ```

Feel free to explore the program's capabilities for ligand-based CTA prediction with your datasets! Please ensure that each dataset is provided as a CSV file, with smiles strings in a column named 'smiles' and compound IDs in a column named 'smiles_id'.

If you use this tool for scientific work that gets published, kindly include a citation of the following paper:

A. Abdulhakeem Mansour Alhasbary and N. Hashimah Ahamed Hassain Malim, "A Similarity-Based Target Prediction Tool for Natural Products," Mol. Inform., vol. -, no. -, p. -, date , doi: .
