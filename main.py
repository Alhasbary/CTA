
"""
Compound-Target Activity (CTA) Prediction Program.

Inputs:
   - List of optional parameters.
   - input: Full path to the data folder containing SMILES string lists (e.g., *_smiles.csv) and the mini-ChEMBL database [Optional]. It may contain more than one list. The tool processes them sequentially. Each file must be a CSV file that contains smiles strings in a column named 'smiles' and compound id in a column named 'smiles_id'.
   - `--output`: Full path to save the results [Optional].
   - output: Full path to save the results [Optional].

Output: 
   - `*_potential_targets_ss.csv`: Contains all likely targets for the chemical compounds in each dataset.
   - `*_CTA_dataset.csv`: Stores the CTA dataset. 

Usage:
   - For help:
     python main.py -h
   - To run the program:
     python main.py [parameters] --input=FullpathToInputFiles --output=FullPathToSaveResults

Usage Example:
    python main.py 
    
Prepared By: Alhasbary
Date: 05/09/2023
"""

import argparse
import glob
import sys
from SharedFunc import *


def parse_args():
    
    desc = f'''Compound-Target Activity (CTA) Prediction Program: 
    - extractDB: Extracting High-Quality Data.
    - processMS: Processing of Molecular Structures.
    - conductSS: Conducting Chemical Similarity Searching.
    - createCTA: Creating Compound-Target Activity Pair Dataset. 
    '''
    print(desc)
    
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("--cs", type=int, default=[8, 9], action="store", dest="confidence_score", nargs='*',
                        help="Confidence score. The confidence scores range from 0, for as yet uncurated data "
                             "entries, to 9, where a single protein target has been assigned a high degree of "
                             "confidence. default=[8, 9]", choices=range(6, 10))
    parser.add_argument("--sv", type=float_range(0.01, 10000.0), default=10000, action="store", dest="standard_value",
                        help="Activity values in nM measurement. default=10000 nM.")
    parser.add_argument("--salt", action="store_false", default=True,
                        help="Whether to disable desalting preprocess. False: disable desalting. default=True")
    parser.add_argument("--charge", action="store_false", default=True,
                        help="Whether to disable neutralization of charge preprocess. False: disable neutralization of "
                             "charge. default=True")
    parser.add_argument("--fingerprint", action="store", default='ecfp', dest="fingerprint",
                        choices=['avalon', 'ecfp', 'fcfp', 'maccs'],
                        help="Desired fingerprint type (avalon, ecfp, fcfp, or maccs). default='ecfp'")
    parser.add_argument("--nBits", action="store", type=int_range(8, 2048), default=2048, dest="nBits",
                        help="Number of bits parameter that specifies the length of the generated (avalon, ecfp, "
                             "or fcfp)fingerprint. default=2048")
    parser.add_argument("--radius", action="store", type=int, default=2, choices=[2, 3], dest="radius",
                        help="Desired radius value for Morgan ECFP/FCFP fingerprints (2, 3). default=2"
                        )
    parser.add_argument("--score", action="store", type=float_range(0.85, 1.0), default=0.85, dest="Tc_score",
                        help="Desired Tc similarity threshold. default=0.85"
                        )
    parser.add_argument("--batch", action="store", type=int_range(16, 512), default=256, dest="batch",
                        help="Desired batch size. To do the searches in chunks based on memory size. default=256"
                        )
    parser.add_argument("--n_jobs", action="store", type=int, default=-1, dest="n_jobs",
                        help="Number of CPU cores to use. default=-1 to use all available CPU cores."
                        )
    parser.add_argument("--agg", action="store", default='median', dest="agg",
                        choices=['min', 'max', 'mean', 'median'],
                        help="Desired aggregation type (min, max, mean, or median). default=median")
    parser.add_argument("--minCompounds", action="store", type=int_range(1, 10000), default=10,
                        dest="minCompounds", help="Number of compounds that specifies the smallest target size in the "
                                                  "compound-target pair data set. default=10")
    parser.add_argument("--maxCompounds", action="store", type=int_range(1, 10000), default=10000,
                        dest="maxCompounds", help="Number of compounds that specifies the largest target size in the "
                                                  "compound-target pair data set. default=10000")
    parser.add_argument("--input", action="store", dest='inputPath', type=str, default='input',
                        help='Full filepath of the data files, which contains the ChEMBL database and the NPs. Default is input.')
    parser.add_argument("--output", action="store", dest='outputPath', type=str,
                        default='output', help='Full filepath, in which the potential targets as well as the CTA data is saved. Default is output')
    args = parser.parse_args()

    return args


def main():
    # get inputs
    args = parse_args()

    # create a folder to save the results if it doesn't exist
    if not os.path.exists(args.outputPath):
        os.makedirs(args.outputPath)
    # begin log
    logPath = os.path.join(args.outputPath, "Preprocessing_log.txt")
    _logFile = open(logPath, 'a')
    log = Tee(sys.stdout, _logFile)

    print("\n", file=log)
    print(" ".join(sys.argv), file=log)
    print("Start time: %s " % (tstamp()), file=log)

    # Retrieved bioactivity data from ChEMBL with predefined selection constraints applied
    selected_ChEMBL = extractDB(log, args)
    
    # Clean SMILES strings, preprocessed by RDKit and MolVS.
    cleaned_ChEMBL = preprocess_chembl(log, args, selected_ChEMBL)

    del selected_ChEMBL
    gc.collect()
    
    valid_chembl_fps, chembl_ids = chembl_fingerprints(log, args, cleaned_ChEMBL)
    
    args.allDB = glob.glob(args.inputPath+"/*smiles.csv")
    if len(args.allDB)==0:
        print("Cannot find SMILES datasets!")
        exit()
    args.allDB = sorted(args.allDB)

    args.datNames = []
    for p in args.allDB:
        args.datNames.append(os.path.basename(p).split("_")[0])

    print(f"Start processing the {len(args.datNames)} datasets: {args.datNames}", file=log)

    for db_name in args.datNames:
        print(f"\n===========\nData set: {db_name}\n===========\n", file=log)
        smiles = os.path.join(args.inputPath, db_name+"_smiles.csv")
        if not os.path.exists(smiles):
            print("Error, file not found", smiles)
            print("It must be a csv file which contains smiles strings in a column named 'smiles' and compound id in "
                  "a column named 'smiles_id'")
            break

        df_db = pd.read_csv(smiles)
        cleaned_smiles = preprocess_smiles(log, args, df_db)
        valid_nps_fps, nps_ids = nps_fingerprints(log, args, cleaned_smiles)
        
        
        # Chemical similarity searches between two data sets.
        ChEMBL_similarity = conductSS(log, args, valid_chembl_fps, chembl_ids , valid_nps_fps, nps_ids)   
        
        if len(ChEMBL_similarity) > 0:
            # Creating potential targets for each smiles string based on similarity searching.
            Potential_targets_path = os.path.join(args.outputPath, f"{db_name}_potential_targets_ss.csv")
            potential_targets_ss(log, args, cleaned_ChEMBL, ChEMBL_similarity, Potential_targets_path)
            
            # Creating compound-target activity pairs for ML training.
            CTA_dataset_path = os.path.join(args.outputPath, f"{db_name}_CTA_dataset.csv")
            createCTA(log, args, cleaned_ChEMBL, ChEMBL_similarity, CTA_dataset_path)
            
    print("\nFinish time: %s " % (tstamp()), file=log)
    _logFile.close()


if __name__ == '__main__':
    main()
