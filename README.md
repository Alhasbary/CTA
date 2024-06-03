# Compound-Target Activity (CTA) Prediction Tool

The Compound-Target Activity (CTA) Prediction Tool is an open-source ligand-based program designed to identify target information associated with chemical compounds. This command-line tool, written in Python and built upon RDKit, employs a two-stage method that utilizes fingerprinting and similarity-based search techniques to identify potential targets.

## Key Features

### 1. **extractDB: Extraction of High-Quality Data from ChEMBL**

This tool allows you to retrieve high-quality bioactivity data from the ChEMBL database. It comes with predefined selection constraints and offers the flexibility to apply user-specified constraints, such as confidence scores and activity values for assays.

### 2. **processMS: Processing of Molecular Structures**

Responsible for cleaning SMILES strings using RDKit and MolVS, this tool performs essential preprocessing steps, including structure normalization, desalting, and charge neutralization. Users have the option to disable desalting and charge neutralization if needed.

### 3. **conductSS: Conducting Chemical Similarity Searching**

Using RDKit, this tool enables chemical similarity searches, a fundamental step in ligand-based CTA prediction. It retrieves potential targets for each given SMILES string. Users can customize fingerprint settings, including type, nBits, and radius parameters that control substructure size. A similarity threshold (Tc) can also be specified.

### 4. **createCTA: Creating Compound-Target Activity Pair Dataset**

The output of the first stage of the tool is a CTA reference dataset, which includes identified targets and their corresponding compounds from ChEMBL. This dataset is valuable for users to analyse and construct models, enabling accurate target predictions for NP compounds.

### 5. **rankTargets: Retrieving ranked potential targets**
The output of the second stage of the tool consists of ranked lists for the chemical compounds in each input query list. These ranked lists contain identified potential targets based on the mean similarity scores of the top k (adjustable parameter option) similar reference compounds. If k is greater than 1, the tool generates potential target lists for each odd number value within the interval [1, k].

## Parameter Settings
The tool provides various parameters to customize the analysis. Below is a detailed description of each parameter, including its name, description, default value, and acceptable range or choices.

| Parameter        | Description                                                                                                                                                    | Default Value  |
|------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------|
| `cs`   | Confidence score (for generate_CTA script). The confidence scores range from 0 (uncurated data) to 9 (high confidence in a single protein target).  | `8`            |
| `sv`    | Activity values in nM measurement (for generate_CTA script).                                                                                       | `10000` nM     |
| `salt`           | Whether to disable desalting preprocess. False: disable desalting.                                                                                            | `True` (enable)         |
| `charge`         | Whether to disable neutralization of charge preprocess. False: disable neutralization of charge.                                                               | `True` (enable)        |
| `fingerprint`    | Desired fingerprint type (avalon, ecfp, fcfp, or maccs).                                                                                                       | `ecfp`         |
| `nBits`          | Number of bits parameter that specifies the length of the generated fingerprint (avalon, ecfp, or fcfp).                                                       | `2048`         |
| `radius`         | Desired radius value for Morgan ECFP/FCFP fingerprints (2 or 3).                                                                                               | `2`  (ECFP4)          |
| `CTA_Tc`         | Desired value for CTA 'Tc' similarity threshold (0.1-1.0).                                                                                                     | `0.85`         |
| `k`          | Desired value for 'top-k' reference compounds (1-11).                                                                                                          | `1`            |
| `batch`          | Desired batch size value (16-512) for chunk-based searches based on memory size.                                                                               | `256`          |
| `n_jobs`         | Number of CPU cores to use.                                                                                                                                   | `-1` (all available CPU cores) |
| `agg`            | Desired aggregation type (min, max, mean, or median).                                                                                                          | `median`       |
| `minCompounds`   | Number of compounds that specifies the smallest target size in the compound-target pair data set.                                                              | `1`            |
| `maxCompounds`   | Number of compounds that specifies the largest target size in the compound-target pair data set.                                                               | `10000`        |
| `inputPath`      | Full filepath of the data files, which contains the ChEMBL database, the NP resource list, and the NP query lists.                                             | `input`        |
| `outputPath`     | Full filepath, in which the potential targets as well as the CTA data is saved.                                                                                | `output`       |

## System Requirements

- Python 3
- General Python Modules: numpy, pandas, joblib
- Special Python Modules: RDKit and MolVS
- ChEMBL SQLite database

## Installation

1. Download the repository and extract files to your desired location.
2. Download the ChEMBL SQLite dataset from [ChEMBL Downloads](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/)
2. Download the COCONUT SMILES dataset from [COCONUT Downloads](https://coconut.naturalproducts.net/download/)
3. Open a command-line window (Windows) or a terminal (Linux).
4. Run the Python scripts as described below.

## Usage

### Create a Custom mini-ChEMBL SQLite Database. 

Create a custom mini-ChEMBL SQLite database tailored to fulfill specific application requirements. In our case, it significantly reduces the storage size from 22.4 GB (ChEMBL32) to just 714 MB.

- Inputs:
   - data: Full path to the ChEMBL dataset [Required].
   - destination: Full path to save the custom ChEMBL version as 'mini_chembl.db' [Optional].

- Output:
   - 'mini_chembl.db' is stored in the specified destination directory, containing customized information for the application's use.

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

### Create a Compound-Target Activity (CTA) dataset

Create a Compound-Target Activity (CTA) dataset to be used as a reference dataset in similarity-based search techniques to identify potential targets. After downloading the COCONUT SMILES dataset (a tab-separated text file), rename the file to NPs_resource.smi.

***Note:*** The steps to create the CTA with default parameter options and datasets take approximately 12 hours to complete. Therefore, this step is performed once using the preferred parameter options and does not need to be rerun unless users wish to change the parameter options or use different datasets/versions than ChEMBL and COCONUT.

- Inputs:
   - List of optional parameters.
   - input: Full path to the data folder containing the NP recource list named "NPs_resource.smi" and the mini-ChEMBL database [Optional]. The NPs_resource file must be a tab-separated text file containing SMILES strings in the first column and compound IDs in the second column.
   - output: Full path to save the results [Optional].

- Output:
   - An updated version of the `mini_chembl.db` file, which includes the CTA dataset.

- Usage:
   - For help:
     ```
     python generate_CTA.py -h
     ```
   - To run the program:
     ```
     python generate_CTA.py [parameters] --input=FullPathToInputFiles --output=FullPathToSaveResults
     ```
  
- Usage Example:
     ```
     python generate_CTA.py
     ```
          
### Identify the potential target(s)

This script takes SMILES-formatted input list(s) of natural product(s) and identifies potential target(s) using fingerprinting and similarity search methods based on the CTA reference dataset.

- Inputs:
   - List of optional parameters.
   - input: Full path to the data folder containing SMILES string lists (e.g., QueryList1_smiles.csv) and the mini-ChEMBL database [Optional]. The folder may contain multiple lists; the tool processes them sequentially. Ensure the file names follow the pattern: QueryListN_smiles.csv, where N is an integer. Each file must be a CSV file with a column named 'smiles' for SMILES strings and a column named 'smiles_id' for compound IDs.
   - output: Full path to save the results [Optional].

- Output: 
   - `QueryListN_potential_targets_based_on_top_k.csv`: Identifies potential targets for the chemical compounds in the query list based on the mean similarity scores of the top k similar reference compounds. If k is greater than 1, the tool generates potential target lists for each odd-numbered value within the interval [1, k].. 
   
- Usage:
   - For help:
     ```
     python predict.py -h
     ```
   - To run the program:
     ```
     python predict.py [parameters] --input=FullpathToInputFiles --output=FullPathToSaveResults
     ```

- Usage Example:
     ```
    python predict.py 
     ```
    
Feel free to explore the program's capabilities for ligand-based CTA prediction with your datasets! Please ensure that each dataset is provided as a CSV file, with smiles strings in a column named 'smiles' and compound IDs in a column named 'smiles_id'.


