"""
This file is part of the Compound-Target Activity (CTA) Prediction Program, used to create the CTA reference dataset.

Inputs:
   - List of optional parameters.
   - input: Full path to the data folder containing the NP recource list named "NPs_resource.csv" and the mini-ChEMBL database [Optional]. The NPs_resource file must be a CSV file containing SMILES strings in a column named 'smiles' and compound IDs in a column named 'smiles_id'.
   - output: Full path to save the results [Optional].

Output:
   - `CTA_dataset.csv`: Stores the CTA dataset.

Usage:
   - For help:
     python generate_CTA.py -h
   - To run the program:
     python generate_CTA.py [parameters] --input=FullPathToInputFiles --output=FullPathToSaveResults

Usage Example:
    python generate_CTA.py

Prepared By: Alhasbary
Date: 05/06/2024
"""

import sys
from SharedFunc import *


def main():
    # get inputs
    args = parse_args()

    # create a folder to save the results if it doesn't exist
    if not os.path.exists(args.outputPath):
        os.makedirs(args.outputPath)
    # begin log
    logPath = os.path.join(args.outputPath, "generate_CTA_log.txt")
    _logFile = open(logPath, 'a')
    log = Tee(sys.stdout, _logFile)

    print("\n", file=log)
    print(" ".join(sys.argv), file=log)
    print("Start time: %s " % (tstamp()), file=log)
    
    if not assert_CTA_resource_files(log, args):
        print("\nFinish time: %s" % (tstamp()), file=log)
        _logFile.close()
        exit()
    
    generate_CTA_reference_dataset(log, args)                  
        
    print("\nFinish time: %s " % (tstamp()), file=log)
    _logFile.close()

if __name__ == '__main__':
    main()

 
