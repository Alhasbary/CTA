
"""
This file is part of the Compound-Target Activity (CTA) Prediction Program, used to create a CTA dataset to be employed as a reference dataset in similarity-based search techniques to identify potential targets. 
After downloading the COCONUT SMILES dataset (a space-separated text file), rename the file to NPs_resource.smi.

Note: The steps to create the CTA with default parameter options and datasets take approximately 12 hours to complete. 
Therefore, this step is performed once using the preferred parameter options and does not need to be rerun unless users wish to change the parameter options or use different datasets/versions than ChEMBL and COCONUT.

Inputs:
   - List of optional parameters.
   - input: Full path to the data folder containing the NP recource list named "NPs_resource.smi" and the mini-ChEMBL database [Optional]. The NPs_resource file must be a space-separated text file containing SMILES strings in the first column and compound IDs in the second column.
   - output: Full path to save the results [Optional].

Output:
   - An updated version of the mini_chembl.db file, which includes the CTA dataset.

Usage:
   - For help:
     python generate_CTA.py -h
   - To run the program:
     python generate_CTA.py [parameters] --input=FullPathToInputFiles --output=FullPathToSaveResults

Usage Example:
    - Download and rename the file COCONUT_DB.smi to NPs_resource.smi and place it in the input folder.
    - Run the program:
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

 
