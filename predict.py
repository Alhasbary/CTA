
"""
Compound-Target Activity (CTA) Prediction Program.
This file is part of the Compound-Target Activity (CTA) Prediction Program, used to identify the potential target(s).

Inputs:
   - List of optional parameters.
   - input: Full path to the data folder containing SMILES string lists (e.g., *_smiles.csv) and the mini-ChEMBL database [Optional]. 
            It may contain more than one list. The tool processes them sequentially. Each file must be a CSV file that contains smiles strings in a column named 'smiles' and compound id in a column named 'smiles_id'.
   - output: Full path to save the results [Optional].

Output: 
   - `*_potential_targets_based_on_top_k.csv`: Contains identify potential targets for the chemical compounds in each dataset based on the mean similarity scores of the top k similar reference compounds. 
                                            If k is greater than 1, the tool generates a potential target list for each odd number value within the interval [1, k]. 
   
Usage:
   - For help:
     python predict.py -h
   - To run the program:
     python predict.py [parameters] --input=FullpathToInputFiles --output=FullPathToSaveResults

Usage Example:
    python predict.py 

Prepared By: Alhasbary
Date: 05/06/2024
"""

import argparse
import sys
from SharedFunc import *


def main():
    # get inputs
    args = parse_args()

    # create a folder to save the results if it doesn't exist
    if not os.path.exists(args.outputPath):
        os.makedirs(args.outputPath)
    # begin log
    logPath = os.path.join(args.outputPath, "predict_log.txt")
    _logFile = open(logPath, 'a')
    log = Tee(sys.stdout, _logFile)

    print("\n", file=log)
    print(" ".join(sys.argv), file=log)
    print("Start time: %s " % (tstamp()), file=log)
    
    if not assert_CTA_resource_files(log, args):
        print("\nFinish time: %s" % (tstamp()), file=log)
        _logFile.close()
        exit()
    
    if not assert_query_resource_files(log, args):
        print("\nFinish time: %s" % (tstamp()), file=log)
        _logFile.close()
        exit()    
            
    # Retrieved the CTA reference dataset 
    desc = f'''\nRetrieving the CTA reference dataset ...'''
    print(desc, file=log)
    
    try:    

        time_s = time.time()
        dbname = args.mini_chembl[0]        
        sqlstr = f'''select DISTINCT
                CTA_smiles.molregno, smiles, uniprot, standard_type, pchembl_value, target_chembl_id
                from
                CTA_smiles, CTA_activity
                where
                CTA_smiles.molregno = CTA_activity.molregno'''                
        CTA = from_db_to_pd(sqlstr, dbname)     
        # Select only necessary columns from CTA dataframe
        columns = ['molregno', 'smiles', 'target_chembl_id', 'uniprot']
        CTA = CTA[columns]                
        # Remove rows that have same info but difference molregno        
        CTA = CTA.drop_duplicates(subset=CTA.columns.difference(['molregno'])).reset_index(drop=True)  
        print(f'Retrieving took {time.time() - time_s:.4f} seconds.', file=log)     
        
        CTA_fps = CTA_gen_fps(log, args, CTA)   

        # Merge 
        CTA = CTA.merge(CTA_fps, on='molregno', how='inner')        
        
    except NameError as e:
        print(f"An error occurred: {e}")
        exit(0)

    print(f"\nStart processing the {args.datNames} Query sets:", file=log)       
    
    for db_name in args.datNames:
        print(f"\n==================\nData set: {db_name}\n==================\n", file=log)
        smiles = os.path.join(args.inputPath, db_name+"_smiles.csv")
        if not os.path.exists(smiles):
            print("Error, file not found", smiles)
            print("It must be a csv file which contains smiles strings in a column named 'smiles' and compound id in "
                  "a column named 'smiles_id'")
            break
    
        # Read the smiles file into a pandas DataFrame
        df_db = pd.read_csv(smiles)
        query_fps = preprocess_Query_SMILES(log, args, df_db)
        
        # Apply the conversion function to the 'fp' column
        query_fps['fp'] = query_fps['fp_base64'].apply(lambda x: convert_from_base64(args, x))          
        
        # Chemical similarity searches between two data sets.
        CTA_similarity = Query_conductSS(log, args, query_fps[['smiles_id', 'fp']], CTA[['molregno', 'target_chembl_id', 'uniprot', 'fp']])      
                    
        if len(CTA_similarity) > 0:            
           
            # Retrieve potential targets for each smiles string according to the mean similarity scores of K nearest neighbors.
            potential_targets_path = os.path.join(args.outputPath, f"using_{args.fingerprint}_the_{db_name}_potential_targets_based_on_top_")
            rankTargets(log, args, CTA_similarity, potential_targets_path)   
            
            del CTA_similarity
            gc.collect()
                     
    print("\nFinish time: %s " % (tstamp()), file=log)
    _logFile.close()


if __name__ == '__main__':
    main()

 