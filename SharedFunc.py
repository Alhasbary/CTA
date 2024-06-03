"""
This file is part of the Compound-Target Activity (CTA) Prediction Program.
    
Prepared By: Alhasbary
Date: 05/06/2024
"""
import argparse
import gc
import glob
from tqdm import tqdm
import numpy as np
import pandas as pd
import csv
import time
import sqlite3
import os

from joblib import Parallel, delayed

from molvs.normalize import Normalizer
from molvs.charge import Reionizer, Uncharger
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem import MACCSkeys, MolFromSmiles, MolToSmiles, DataStructs
from rdkit.Avalon.pyAvalonTools import GetAvalonFP
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from rdkit.DataStructs.cDataStructs import ExplicitBitVect

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


def parse_args():           
    parser = argparse.ArgumentParser()
    parser.add_argument("--cs", type=int, default=8, action="store", dest="confidence_score", nargs='*',
                        help="Confidence score (for generate_CTA script). The confidence scores range from 0, for as yet uncurated data "
                             "entries, to 9, where a single protein target has been assigned a high degree of "
                             "confidence. default=8", choices=range(6, 10))
    parser.add_argument("--sv", type=float_range(0.01, 10000.0), default=10000, action="store", dest="standard_value",
                        help="Activity values in nM measurement (for generate_CTA script). default=10000 nM.")
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
                        help="Desired radius value for Morgan ECFP/FCFP fingerprints (2, 3). default=2")
    parser.add_argument("--CTA_Tc", action="store", type=float_range(0.1, 1.0), default=0.85, dest="CTA_Tc",
                        help="Desired value for CTA 'Tc' similarity threshold (0.1-1.0). default=0.85.")
    parser.add_argument("--k", action="store", type=int_range(1, 11), default=3, dest="top_k",
                        help="Desired value for 'top-k' reference compounds (1-11). Default is 3.")
    parser.add_argument("--batch", action="store", type=int_range(16, 512), default=256, dest="batch",
                        help="Desired batch size value (16-512) for chunk-based searches based on memory size. Default is 256.")
    parser.add_argument("--n_jobs", action="store", type=int, default=-1, dest="n_jobs",
                        help="Number of CPU cores to use. default=-1 to use all available CPU cores.")
    parser.add_argument("--agg", action="store", default='median', dest="agg",
                        choices=['min', 'max', 'mean', 'median'],
                        help="Desired aggregation type (min, max, mean, or median). default=median")
    parser.add_argument("--minCompounds", action="store", type=int_range(1, 10000), default=1,
                        dest="minCompounds", help="Number of compounds that specifies the smallest target size in the "
                                                  "compound-target pair data set. default=1")
    parser.add_argument("--maxCompounds", action="store", type=int_range(1, 10000), default=10000,
                        dest="maxCompounds", help="Number of compounds that specifies the largest target size in the "
                                                  "compound-target pair data set. default=10000")
    parser.add_argument("--input", action="store", dest='inputPath', type=str, default='input',
                        help='Full filepath of the data files, which contains the ChEMBL database, the NP recource list, and the NP query lists. Default is input.')
    parser.add_argument("--output", action="store", dest='outputPath', type=str,
                        default='output', help='Full filepath, in which the potential targets as well as the CTA data is saved. Default is output')
    args = parser.parse_args()
    return args
    
def assert_CTA_resource_files(log, args):
    """
    Make sure the existence of the reseource files
    """
    args.mini_chembl = glob.glob(args.inputPath + "/mini_chembl.db")
    if len(args.mini_chembl) == 0:
        print("Cannot find the mini-ChEMBL database file!", file=log)
        print("The mini-ChEMBL database file must be an SQLite database named 'mini_chembl.db'. Please refer to 'README.md' for more information.", file=log)            
        return False

    args.NPs_resource = glob.glob(args.inputPath + "/NPs_resource.smi")
    if len(args.NPs_resource) == 0:
        print("Cannot find the NPs resource file!", file=log)
        print("The NPs resource file must be a CSV named 'NPs_resource.csv' containing SMILES strings in a column named 'smiles' and compound IDs in a column named 'smiles_id'.", file=log)            
        return False
        
    return True
                
def assert_query_resource_files(log, args):
    """
    Make sure the existence of the query list files
    """
    args.mini_chembl = glob.glob(args.inputPath + "/mini_chembl.db")
    if len(args.mini_chembl) == 0:
        print("Cannot find the mini-ChEMBL database file!", file=log)
        print("The mini-ChEMBL database file must be an SQLite database named 'mini_chembl.db'. Please refer to 'README.md' for more information.", file=log)            
        return False

    args.allQueryLists = glob.glob(args.inputPath+"/*smiles.csv")
    if len(args.allQueryLists)==0:
        print("Cannot find any query list file!", file=log)
        print("The query list file must be a CSV file named '*_smiles.csv' containing SMILES strings in a column named 'smiles' and compound id in a column named 'smiles_id'. The tool ccan processes more than one query list.", file=log)            
        return False
       
    args.allQueryLists = sorted(args.allQueryLists)
    args.datNames = []
    for p in args.allQueryLists:
        args.datNames.append(os.path.basename(p).split("_")[0])
    return True
                    
def tstamp():
    t = time.localtime()
    return "%d/%d/%d %d:%.2d.%d" % ((t[1], t[2], t[0]) + t[3:6])

class Tee(object):
    def __init__(self, *args):
        self.files = [f for f in args]

    def write(self, s):
        for f in self.files:
            f.write(s)

    def flush(self):
        for f in self.files:
            f.flush()

def print_statistics(log, df):
    print('Compound-target pairs:', file=log)
    print(f"    Number of bioactivity records    = {len(df)}", file=log)
    print(f"    Number of unique uniprot targets = {len(df.uniprot.unique())}", file=log)
    print(f"    Number of unique compounds       = {len(df.smiles.unique())}", file=log)
    #print(f"    Target types included: {df.target_type.unique()}", file=log)
    df_a = df.groupby(['uniprot']).count()
    print('\nCompound statistics in the targets:', file=log)
    print('min = ', df_a.smiles.min(), '\t median = ', df_a.smiles.median(), '\t mean = ', df_a.smiles.mean(),
          '\t max = ', df_a.smiles.max(), '\n', file=log)

def int_range(min_value, max_value):
    def int_range_type(value):
        ivalue = int(value)
        if ivalue < min_value or ivalue > max_value:
            raise argparse.ArgumentTypeError("Value must be between %s and %s. For the Avalon fingerprint, the value "
                                             "is recommended to be divisible by 8." % (min_value, max_value))
        return ivalue

    return int_range_type

def float_range(min_value, max_value):
    def float_range_type(value):
        fvalue = float(value)
        if fvalue < min_value or fvalue > max_value:
            raise argparse.ArgumentTypeError("Value must be between %s and %s" % (min_value, max_value))
        return fvalue
    return float_range_type


def from_db_to_pd(sql, dbname):
    """
    Retrieved bioactivity data from ChEMBL SQLite database into a DataFrame object    
    :param sql: String, the SQL statement.
    :param dbname: String, Database name
    :return: Dataframe contains the selected bioactivity data.
    """
    conn = sqlite3.connect(dbname)
    df = pd.read_sql(sql, conn)
    return df
  
def extractDB(log, args):
    """
    Extract High-Quality Data.
    :param log: Log file for recording extraction information.
    :param args: argparse object containing program arguments.
    :return: DataFrame containing the selected data.
    """
    try:
        sqlstr = f'''select DISTINCT
            cs.molregno, md.chembl_id as molecule_chembl_id, cs.canonical_smiles, act.standard_value, 
            act.standard_type, act.pchembl_value, ass.confidence_score, ass.description,
            td.organism, td.target_type, td.pref_name, td.chembl_id as target_chembl_id, 
            cseq.accession as uniprot, pc.protein_class_desc, cp.np_likeness_score 
            from
            molecule_dictionary md, compound_structures cs, 
            activities act, assays ass, target_dictionary td, target_components tc,
            component_sequences cseq, component_class cclass, protein_classification pc, compound_properties cp
            where
            cs.molregno = md.molregno and
            cp.molregno = md.molregno and
            md.molregno = act.molregno and
            act.assay_id = ass.assay_id and
            ass.tid = td.tid and 
            tc.tid = td.tid and
            cseq.component_id = tc.component_id and
            cseq.component_id = cclass.component_id and
            cclass.protein_class_id = pc.protein_class_id and
            ass.confidence_score >= {args.confidence_score} and
            act.standard_value <= {args.standard_value}'''        

        dbname = os.path.join(args.inputPath, "mini_chembl.db")
        selected_bioactivity = from_db_to_pd(sqlstr, dbname)
        selected_bioactivity = selected_bioactivity.drop_duplicates(
            subset=selected_bioactivity.columns.difference(['protein_class_desc']))
        return selected_bioactivity

    except NameError as e:
        print(f"An error occurred: {e}")
        exit(0)
        
def processMS(smi_id, smi, args):
    """
    Process Molecular Structures.
    :param smi: SMILES string.
    :param args: argparse object containing program arguments.
    :return: List of SMILES strings with standardized structures. Problematic SMILES strings are replaced with 'nan' and removed.
    """  
    fp = 'issue'    
    normalizer = Normalizer()
    remover = SaltRemover()
    neutralize1 = Reionizer()
    neutralize2 = Uncharger()
    mol = MolFromSmiles(smi)
    if mol is not None:
        mol = normalizer.normalize(mol)
        if args.salt:
            mol = remover(mol)
        if args.charge:
            mol = neutralize1(mol)
            mol = neutralize2(mol)
        smiles = MolToSmiles(mol)
        try:
            if args.fingerprint == 'avalon':
                fp = GetAvalonFP(mol, nBits=args.nBits)
            elif args.fingerprint == 'ecfp':
                fp = GetMorganFingerprintAsBitVect(mol, nBits=args.nBits, radius=args.radius)
            elif args.fingerprint == 'fcfp':
                fp = GetMorganFingerprintAsBitVect(mol, nBits=args.nBits, radius=args.radius, useFeatures=True)
            elif args.fingerprint == 'maccs':
                fp = MACCSkeys.GenMACCSKeys(mol)
        except:
            # print(f'Issue with conversion to {args.fingerprint} fingerprint: ' + str(smi))
            pass
    else:
        smiles = 'issue'
    if fp != 'issue':
        fp = fp.ToBase64()
    return smi_id, smiles, fp

def gen_fp(smi, args):
    """
    Generate a fingerprint for a given SMILES string based on specified fingerprint type.    
    Parameters:
    smi (str): SMILES string representing the molecule.
    args (Namespace): Arguments containing fingerprint type and parameters.    
    Returns:
    str: Base64 encoded fingerprint or 'issue' if an error occurs.
    """
    # Convert SMILES to RDKit molecule
    mol = MolFromSmiles(smi)
    if mol is None:
        return 'issue'        
    try:
        # Generate the appropriate fingerprint
        if args.fingerprint == 'avalon':
            fp = GetAvalonFP(mol, nBits=args.nBits)
        elif args.fingerprint == 'ecfp':
            fp = GetMorganFingerprintAsBitVect(mol, nBits=args.nBits, radius=args.radius)
        elif args.fingerprint == 'fcfp':
            fp = GetMorganFingerprintAsBitVect(mol, nBits=args.nBits, radius=args.radius, useFeatures=True)
        elif args.fingerprint == 'maccs':
            fp = MACCSkeys.GenMACCSKeys(mol)
        else:
            return 'issue'  # Handle unknown fingerprint type
    except Exception as e:
        # Log the error or handle it as needed
        # print(f'Issue with conversion to {args.fingerprint} fingerprint: {str(smi)} - {str(e)}')
        return 'issue'    
    # Encode the fingerprint in Base64
    return fp.ToBase64()    
        
def CTA_gen_fps(log, args, CTA):
    CTA = CTA.drop_duplicates(subset=['smiles']).reset_index(drop=True) 
    print(f'Generating fingerprints for {len(CTA)} CTA smiles ...', file=log)
    time_s = time.time()
    fps = Parallel(n_jobs=int(args.n_jobs), backend="multiprocessing")(
        delayed(gen_fp)(row['smiles'], args) for _, row in tqdm(CTA.iterrows(), total=len(CTA), desc=" gen_fp"))  
    
    CTA_fps = pd.DataFrame()
    CTA_fps['fp_base64'] = fps
    CTA_fps['molregno'] = CTA['molregno']        
    # Remove rows containing problematic compounds
    CTA_fps = CTA_fps[CTA_fps['fp_base64'] != 'issue'].reset_index(drop=True)     
    CTA_fps['fp'] = CTA_fps['fp_base64'].apply(lambda x: convert_from_base64(args, x))
    
    print(f'Generating fingerprints took {time.time() - time_s: .4f} seconds.',
          file=log)
    return CTA_fps[['molregno', 'fp']]
        
def generate_CTA_reference_dataset(log, args):
    """
    Generate the CTA reference dataset
    :param log: Log file for recording processing information.
    :param args: argparse object containing program arguments.
    """
    try:
        dbname = args.mini_chembl[0]       
        mini_chembl = sqlite3.connect(dbname)
        mini_chembl_cur = mini_chembl.cursor()
        
        sqlstr = f'''select DISTINCT confidence_score, standard_value, salt, charge, fingerprint , nBits, radius, CTA_Tc from init'''    
        parameters = from_db_to_pd(sqlstr, dbname)
        confidence_score = parameters['confidence_score'].values[0]
        standard_value = parameters['standard_value'].values[0]
        salt = parameters['salt'].values[0]
        charge = parameters['charge'].values[0]
        fingerprint = parameters['fingerprint'].values[0]
        nBits = parameters['nBits'].values[0]
        radius = parameters['radius'].values[0]
        CTA_Tc = parameters['CTA_Tc'].values[0]
        desc = f'''Compound-Target Activity (CTA) Prediction Program:
                Steps of Generating the CTA:
                1. Data Extraction from ChEMBL:
                    Selects bioactivity data based on specific constraints:
                        - Canonical SMILES (not null)
                        - Accession (not null)
                        - Assay type (binding assays)
                        - Confidence score >= {confidence_score} (Parameter option: --cs)
                        - Standard units (nM)
                        - Standard value <= {standard_value} (Parameter option: --sv)
                        - Standard relation (=)
                        - Standard type (IC50, Ki, Kd, or EC50)
                        - pChEMBL value (not null)
                        - Data validity comment (null or manually validated)
                        - Activity comment (not undetermined, inconclusive, unspecified, inactive, or not active)
                        - Potential duplicates (0)

                2. ChEMBL and NPs SMILES Preprocessing:
                    Utilizes RDKit and MolVS to preprocess SMILES strings, including:
                        - Desalting (Parameter option: --salt)
                        - Charge neutralization (Parameter option: --charge)
                    Users can disable desalting and charge neutralization using the respective options.

                3. Generation of Molecular Fingerprints:
                    Employs RDKit to create molecular fingerprints for each processed SMILES string, considering:
                        - Fingerprint type ({fingerprint}): Options include avalon, ecfp, fcfp, or maccs (Parameter option: --fingerprint)
                        - Fingerprint length ({nBits}) (Parameter option: --nBits)
                        - Fingerprint radius ({radius}) (Parameter option: --radius)

                4. Generate CTA Data Set:
                    Performs chemical similarity searches, considering:
                        - ChEMBL records whose targets interact with at least one compound having a Tc of {CTA_Tc} or higher with any of the source NP compounds are considered and saved (Parameter option: --CTA_Tc)

                To adjust the parameter options for the generating steps, rerun the program and specify the desired options. 
                The processing time for the generating steps with default parameter options is approximately 12 hours.

                For a list of available parameter options, run the program with the -h argument.
                Here are the default values for some other parameter options:
                    - Aggregation type ({args.agg}): Options include min, max, mean, or median (Parameter option: --agg)
                    - Minimum number of compounds ({args.minCompounds}) (Parameter option: --minCompounds)
                    - Maximum number of compounds ({args.maxCompounds}) (Parameter option: --maxCompounds)\n
                '''
        if CTA_Tc == 100:
            initial_run = True
        else:
            initial_run = False
            
        if not initial_run:    
            if (confidence_score!= args.confidence_score) or \
            (standard_value != args.standard_value) or \
            (salt != int(args.salt)) or \
            (charge != int(args.charge)) or \
            (fingerprint != args.fingerprint) or \
            (nBits != args.nBits) or \
            (radius != args.radius) or \
            (CTA_Tc != args.CTA_Tc):        
                print(desc, file=log)      
                print('\nThe steps to creat the CTA with default parameter options take approximately 12 hours to complete.')
                inp = input('Would you like to continue with the new parameter options? ([y]/n): ').strip().lower() or 'y'
              
            else:
                print('\nThe CTA reference dataset has already been created with the given parameter options!', file=log)                  
                inp = 'n'            
            
        if initial_run or inp == 'y':  
            print('\nStarting the generating steps. Estimated completion time is around 12 hours.', file=log)      
            
            # Preprocess all ChEMBL SMILES strings, clean and generate fingerprints by RDKit and MolVS.
            clean_ChEMBL_structures = preprocess_ChEMBL_SMILES(log, args)  
            ChEMBL_fps = clean_ChEMBL_structures[['uniprot', 'molregno', 'fp_base64']].drop_duplicates().reset_index(drop=True)
            ChEMBL_fps['fp'] = ChEMBL_fps['fp_base64'].apply(lambda x: convert_from_base64(args, x))      
                    
            # Preprocess all NPs SMILES strings, clean and generate fingerprints by RDKit and MolVS.
            fname = args.NPs_resource[0]                
            clean_NPs_structures = preprocess_NPs_SMILES(log, args, fname)        
            clean_NPs_structures['fp'] = clean_NPs_structures['fp_base64'].apply(lambda x: convert_from_base64(args, x))    
            
            # Extracting unique UniProt identifiers for ChEMBL targets interacting with similar compounds
            interact_ChEMBL_targets = CTA_conductSS(log, args, clean_NPs_structures[['smiles_id', 'fp']], ChEMBL_fps[['uniprot', 'molregno', 'fp']], args.CTA_Tc)  
            
            if len(interact_ChEMBL_targets) > 0:                                                    
            
                # Creating a subset of ChEMBL structures data containing only interactions with identified targets
                partial_ChEMBL = clean_ChEMBL_structures[clean_ChEMBL_structures['uniprot'].isin(interact_ChEMBL_targets)].reset_index(drop=True)
                
                # Save the partial_ChEMBL to a CSV file
                #partial_ChEMBL.to_csv(os.path.join(args.outputPath,'subset_of_ChEMBL.csv'), index=False, float_format='%.4f')
                
                # Creating compound-target activity pairs.      
                CTA = createCTA(log, args, partial_ChEMBL)      
                    
                # Update CTA data set table                    
                mini_chembl_cur.executescript('''
                DROP TABLE IF EXISTS CTA_smiles;
            
                CREATE TABLE CTA_smiles (
                molregno BIGINT NOT NULL,
                smiles VARCHAR(4000),
                CONSTRAINT pk_molregno_CTA PRIMARY KEY (molregno)
                );''')                            
                
                CTA_smiles = CTA.drop_duplicates(subset=['molregno']).reset_index(drop=True)
                for idx, row in CTA_smiles.iterrows():
                    mini_chembl_cur.execute('''INSERT INTO CTA_smiles (molregno, smiles)
                     VALUES (?, ?)''', (row['molregno'], row['smiles']))
                
                mini_chembl_cur.executescript('''
                DROP TABLE IF EXISTS CTA_activity;
            
                CREATE TABLE CTA_activity (
                molregno BIGINT NOT NULL,
                uniprot VARCHAR(25),
                standard_type VARCHAR(250),
                pchembl_value NUMERIC(4, 2),
                target_chembl_id VARCHAR(20) NOT NULL,
                CONSTRAINT fk_molregno_CTA FOREIGN KEY(molregno) REFERENCES CTA_smiles (molregno) ON DELETE CASCADE
                );''')        
                                    
                for idx, row in CTA.iterrows():
                    mini_chembl_cur.execute('''INSERT INTO CTA_activity (molregno, uniprot, standard_type, pchembl_value, target_chembl_id)
                     VALUES (?, ?, ?, ?, ?)''', (row['molregno'], row['uniprot'], row['standard_type'], row['pchembl_value'], row['target_chembl_id']))
                     
                # Update parameter setting table                 
                mini_chembl_cur.execute('''UPDATE  init SET  confidence_score=?, standard_value=?, 
                salt=?, charge=?, fingerprint=? , nBits=?, radius=?, CTA_Tc=?''', (args.confidence_score,args.standard_value, int(args.salt),int(args.charge),
                                                                                   args.fingerprint,args.nBits,args.radius, args.CTA_Tc)) 
                                                     
                mini_chembl.commit()
                mini_chembl_cur.close()
                mini_chembl.close()
            
    
    except NameError as e:
        print(f"An error occurred: {e}")
        exit(0)
        
def preprocess_ChEMBL_SMILES(log, args):
    """
    Preprocess ChEMBL Molecular Structures.
    :param log: Log file for recording processing information.
    :param args: argparse object containing program arguments.
    """
    try:
        desc = f'''\nExtract and parallel preprocess ChEMBL SMILES strings: clean and generate fingerprints by RDKit and MolVS.'''
        print(desc, file=log)
        time_s = time.time()
                
        # Retrieved bioactivity data from ChEMBL with predefined selection constraints applied        
        ChEMBL_info = extractDB(log, args)
        # Save related information
        columns = ['molregno', 'canonical_smiles', 'standard_type', 'pchembl_value', 'target_chembl_id', 'target_type', 'uniprot']   
        ChEMBL_info = ChEMBL_info[columns]#[:20000]    
        ChEMBL_info.rename(columns={'canonical_smiles': 'smiles'}, inplace=True)        
        print('ChEMBL info after applying the selection protocol:', file=log)
        print_statistics(log, ChEMBL_info)
             
        ChEMBL_unique = ChEMBL_info.drop_duplicates(subset=['smiles']).reset_index(drop=True)       
        cleaned_smi = Parallel(n_jobs=int(args.n_jobs), backend="multiprocessing", batch_size=args.batch)(
            delayed(processMS)(row['molregno'], row['smiles'], args) for _, row in tqdm(ChEMBL_unique.iterrows(), total=len(ChEMBL_unique), desc="ChEMBL processMS")) 
            
        clean_structures = pd.DataFrame(cleaned_smi, columns=['molregno', 'smiles', 'fp_base64'])            
        ChEMBL_info.drop('smiles', axis=1, inplace=True)
        clean_ChEMBL_structures = pd.merge(ChEMBL_info, clean_structures, on='molregno', how='inner')
        
        # Remove rows containing problematic compounds
        clean_ChEMBL_structures = clean_ChEMBL_structures[clean_ChEMBL_structures['fp_base64'] != 'issue'].reset_index(drop=True)       
        print('ChEMBL info after preprocessing:', file=log)
        print_statistics(log, clean_ChEMBL_structures)
        
        elapsed_time = time.time() - time_s
        print(f'Preprocess ChEMBL SMILES strings completed in {elapsed_time:.2f} seconds.', file=log)          
        return clean_ChEMBL_structures
        
    except NameError as e:
        print(f"An error occurred: {e}")
        exit(0)
        
def preprocess_NPs_SMILES(log, args, fname):
    """
    Preprocess NPs Molecular Structures.
    :param log: Log file for recording processing information.
    :param args: argparse object containing program arguments.
    """
    try:
        desc = f'''\nExtract and parallel preprocess the NP scource compounds: clean and generate fingerprints by RDKit and MolVS.'''
        print(desc, file=log)
        time_s = time.time()        
        
        NPs_info = pd.read_csv(fname, sep=' ', header=None, names=['smiles', 'smiles_id']) 
        NPs_info = NPs_info.drop_duplicates(subset=['smiles']).reset_index(drop=True)
        
        print('Number of compounds in the NP scource file:', NPs_info.shape[0], file=log)
        NPs_info= NPs_info#[:10000]
        
        cleaned_smi = Parallel(n_jobs=int(args.n_jobs), backend="multiprocessing", batch_size=args.batch)(
            delayed(processMS)(row['smiles_id'], row['smiles'], args) for _, row in tqdm(NPs_info.iterrows(), total=len(NPs_info), desc="NPs processMS"))  
            
        clean_structures = pd.DataFrame(cleaned_smi, columns=['smiles_id', 'smiles', 'fp_base64'])    
        NPs_info.drop('smiles', axis=1, inplace=True)
        clean_NPs_structures = pd.merge(NPs_info, clean_structures, on='smiles_id', how='inner')
                
        # Remove rows containing problematic compounds
        clean_NPs_structures = clean_NPs_structures[clean_NPs_structures['fp_base64'] != 'issue'].reset_index(drop=True) 
        print('clean NPs structures:', clean_NPs_structures.shape[0], file=log)
        
        elapsed_time = time.time() - time_s
        print(f'Preprocess NPs SMILES strings completed in {elapsed_time:.2f} seconds.', file=log)
        
        return clean_NPs_structures
        
    except NameError as e:
        print(f"An error occurred: {e}")
        exit(0)
    
def preprocess_Query_SMILES(log, args, df_db):
    """
    Process NPs Molecular Structures.
    :param log: Log file for recording processing information.
    :param args: argparse object containing program arguments.
    :param df_db: DataFrame containing the NPs data.
    :return: DataFrame containing the cleaned NPs data.
    """
    try:
        print('SMILES strings before preprocessing:', file=log)
        print(f'Number of compounds: {len(df_db.smiles_id)}', file=log)
        df_db = df_db.drop_duplicates(subset=['smiles']).reset_index(drop=True) 
        print(f'Number of unique compounds: {len(df_db.smiles_id.unique())}', file=log)

        # Start clean preprocess in parallel
        time_s = time.time()
        cleaned_smi = Parallel(n_jobs=int(args.n_jobs), backend="multiprocessing", batch_size=args.batch)(
            delayed(processMS)(row['smiles_id'], row['smiles'], args) for _, row in tqdm(df_db.iterrows(), total=len(df_db), desc="Query processMS"))               
                
        clean_Query_structures = pd.DataFrame(cleaned_smi, columns=['smiles_id', 'smiles', 'fp_base64'])
                        
        # Identify and handle rows with problematic compounds
        smiles_id_issue = clean_Query_structures[clean_Query_structures['fp_base64'] == 'issue']['smiles_id'].values
        if len(smiles_id_issue) > 0:
            print(f'Unable to generate {args.fingerprint} fingerprint for {len(smiles_id_issue)} SMILES strings:', file=log)
            print(list(smiles_id_issue), file=log)
        cleaned_data = clean_Query_structures[clean_Query_structures['fp_base64'] != 'issue'].reset_index(drop=True)

        # Print statistics and save the results
        print('\nSMILES strings after preprocessing:', file=log)
        print(f'Number of compounds: {len(cleaned_data.smiles_id)}', file=log)
        print(f'Number of unique compounds: {len(cleaned_data.smiles_id.unique())}', file=log)
        elapsed_time = time.time() - time_s
        print(f' preprocessing completed in {elapsed_time:.2f} seconds.', file=log)
        return cleaned_data
        
    except NameError as e:
        print(f"An error occurred: {e}")
        exit(0)

def convert_from_base64(args, x):
    fp_from_base64 = ExplicitBitVect(args.nBits)
    fp_from_base64.FromBase64(x)
    return fp_from_base64
    
def ss(args, row, reference_fps, Tc):
    ss_results = pd.Series(DataStructs.BulkTanimotoSimilarity(row['fp'], list(reference_fps['fp'])), index=reference_fps.index)
    mask = ss_results >= Tc
    return reference_fps[mask]['uniprot'].unique()
    
def CTA_conductSS(log, args, query_fps, reference_fps, Tc):
    desc = f'''
    Perform chemical similarity searches. 
    The cleaned SMILES strings for each compound are converted into {args.nBits}-bit {args.fingerprint} fingerprints 
    with a radius of {args.radius}. From {len(reference_fps.uniprot.unique())} targest, targets that interact with at least one compound having a Tc of {args.CTA_Tc} 
    or higher with any of the scource compounds are considered and saved.
    '''
    
    print(desc, file=log)
    query_fps = query_fps.drop_duplicates()    
    time_s = time.time()
       
    # Showing the progress bar increases the estimated completion time to around 17 hours.
    '''
    results = Parallel(n_jobs=args.n_jobs, backend="multiprocessing", batch_size=args.batch)(
        delayed(ss)(args, row, reference_fps, Tc) for _, row in tqdm(query_fps.iterrows(), total=len(query_fps), desc="CTA conductSS"))    
    '''

    # Estimated completion time without the progress bar is around 12 hours.
    print('Estimated completion time without the progress bar is around 12 hours.')
    print('Showing the progress bar increases the estimated completion time to around 17 hours.\n')
    results = Parallel(n_jobs=args.n_jobs, backend="multiprocessing", batch_size=args.batch)(
        delayed(ss)(args, row, reference_fps, Tc) for _, row in query_fps.iterrows())
        
        
    # Flatten the list of lists and get unique UniProt IDs
    uniprot_ids = set(target for sublist in results for target in sublist)
    
    elapsed_time = time.time() - time_s
    if uniprot_ids:
        print(f'Similarity searches completed in {elapsed_time:.2f} seconds. The mini_chembl.db file has been updated with the generated CTA dataset.', file=log)
    else:
        print(f'No similarities were found!. The process completed in {elapsed_time:.2f} seconds.', file=log)

    return list(uniprot_ids)
    
def Query_conductSS(log, args, query_fps, reference_fps):
    desc = f'''\nPerform chemical similarity searches. 
    The cleaned SMILES strings for each compound are converted into {args.nBits}-bit {args.fingerprint} fingerprints 
    with a radius of {args.radius}.'''    
    
    print(desc, file=log)
    query_fps = query_fps.drop_duplicates()
    
    print(f'\nStart chemical similarity searches for {len(query_fps)} smiles\n', file=log)
    time_s = time.time()

    # Initialize lists to store relevant information
    smiles_ids = []
    molregnos = []
    target_chembl_ids = []
    uniprots = []
    scores = []
    
    # Group reference_fps by 'uniprot'
    grouped_fps = reference_fps.groupby('uniprot')
    
    # Iterate over each query
    for idx, row in tqdm(query_fps.iterrows(), total=len(query_fps), desc="Query conductSS"):
        
        # Calculate similarity scores for each group (target)
        for target, fps in grouped_fps:
            # Calculate similarity scores
            ss_results = pd.Series(DataStructs.BulkTanimotoSimilarity(row['fp'], list(fps['fp'])), index=fps.index)
            
            # Get indices of top N similarity scores
            top_indices = ss_results.nlargest(args.top_k, keep='all').index
            
            # Store relevant information
            smiles_ids.extend([row['smiles_id']] * len(top_indices))
            molregnos.extend(fps.loc[top_indices, 'molregno'].tolist())
            target_chembl_ids.extend(fps.loc[top_indices, 'target_chembl_id'].tolist())
            uniprots.extend(fps.loc[top_indices, 'uniprot'].tolist())
            scores.extend(ss_results[top_indices].tolist())

    # Create DataFrame with relevant information
    df_result = pd.DataFrame({
        'smiles_id': smiles_ids,
        'molregno': molregnos,
        'target_chembl_id': target_chembl_ids,
        'uniprot': uniprots,
        'score': scores
    })

    if not df_result.empty:
        elapsed_time = time.time() - time_s
        print(f'Similarity search completed in {elapsed_time:.2f} seconds.', file=log)     
    else:
        print('No similarities were found!', file=log)

    return df_result


def createCTA(log, args, selected_ChEMBL):
    """
    Creating a Compound-Target Activity (CTA) Pair Dataset.

    :param log: Log file for recording dataset creation information.
    :param args: argparse object containing program arguments.
    :param selected_ChEMBL: DataFrame containing the ChEMBL dataset with target information.
    :param dest: File name for saving the CTA dataset.
    """
    
    desc = f'''\nCreate a compound-target pair dataset. Only compounds with a similarity score of {args.CTA_Tc} or 
    higher are considered from the cleaned ChEMBL dataset. When identical target-compound pairs are found, the {args.agg} 
    activity value is used. Additionally, the dataset is pruned, retaining only compound-target pairs where targets have 
    at least {args.minCompounds} compounds. Furthermore, the number of compounds per target is capped at {args.maxCompounds} through random sampling.\n'''
    print(desc, file=log)
    time_s = time.time()
    
    if selected_ChEMBL.empty:
        print(f'No compound having a Tanimoto coefficient (Tc) of {args.CTA_Tc} or higher were found!', file=log)
        return []
    
    # When identical target-compound pairs are found, the {args.agg} activity value is used.
    column_map = {col: "first" for col in selected_ChEMBL.columns}
    column_map["pchembl_value"] = args.agg
    CTA_multiple_activity = selected_ChEMBL.groupby(["uniprot", "standard_type", "smiles"],
                                             as_index=False).aggregate(column_map)        
    
    # Retain only labels with at least args.minCompounds observations
    subset_lower = CTA_multiple_activity.groupby('uniprot').uniprot.transform(len) >= args.minCompounds
    CTA_pruned = CTA_multiple_activity.loc[subset_lower]

    # Find labels with more than args.maxCompounds examples each
    subset_gt = CTA_pruned.groupby('uniprot').uniprot.transform(len) > args.maxCompounds
    CTA_gt = CTA_pruned.loc[subset_gt]
    CTA_lt = CTA_pruned.loc[~subset_gt]
    
    # Keep only args.maxCompounds examples of those labels
    subset_upper = CTA_gt.groupby('uniprot', group_keys=True).apply(pd.DataFrame.sample,
                                                                     n=args.maxCompounds).reset_index(level='uniprot',
                                                                                                      drop=True)    
    
    CTA = pd.concat([CTA_lt, subset_upper], axis=0)
    print('\nCTA data set:', file=log)
    print_statistics(log, CTA)    
    
    #CTA.sort_values(by=['molregno', 'pchembl_value'], inplace=True, ascending=[True, False])    
    #CTA.to_csv(os.path.join(args.outputPath, "CTA_dataset.csv"), index=False, float_format='%.4f')
    
    print(f'\nCreating the compound-target activity pair dataset took {time.time() - time_s:.4f} seconds.', file=log)
    return CTA

def rankTargets(log, args, reference_dataset, dest):
    """
    Retrieving ranked potential targets.
    :param log: Log file for recording search information.
    :param args: argparse object containing program arguments.
    :param reference_dataset: DataFrame containing the results of conductSS.
    :param dest: File name for saving the potential targets for each SMILES string.
    """ 
          
    desc = f'\nRetrieve potential targets for each SMILES string.\n'
    print(desc, file=log)
    time_s = time.time()    
    
    if reference_dataset.empty:
        print(f'The provided CTA reference dataset is empty!', file=log)
        return 
    
    # List of k-values to process
    if args.top_k > 1:
        k_values = [i for i in range(1, args.top_k + 1) if i % 2 != 0]
    else:
        k_values = [1]
        
    # Group reference_dataset by ['uniprot', 'smiles_id']
    grouped = reference_dataset.groupby(['uniprot', 'smiles_id'])

    # Loop over each k-value
    for top_k in k_values:
        # Initialize lists to store relevant information for each k-value
        smiles_ids = []
        target_chembl_ids = []
        uniprots = []
        mean_scores = []

        # Iterate over each group (['uniprot', 'smiles_id'])
        for _, group_df in tqdm(grouped, desc=f"k = {top_k}", total=len(grouped)):
            # Select the top K scores within the group
            top_K_group = group_df.nlargest(top_k, 'score', keep='all')

            # Calculate the mean of the top K scores for the group
            mean_score = top_K_group['score'].mean()

            # Append relevant information to lists
            smiles_ids.append(group_df['smiles_id'].iloc[0])
            target_chembl_ids.append(group_df['target_chembl_id'].iloc[0])
            uniprots.append(group_df['uniprot'].iloc[0])
            mean_scores.append(mean_score)

        # Create DataFrame with relevant information for the current k-value
        df = pd.DataFrame({
            'smiles_id': smiles_ids,
            'target_chembl_id': target_chembl_ids,
            'uniprot': uniprots,
            'mean_score': mean_scores,
        })

        # Calculate ranks based on the mean score, handling each 'smiles_id' separately
        df['rank'] = df.groupby('smiles_id')['mean_score'].rank(ascending=False, method='min').astype(int)

        # Sort values by smiles_id and mean_score
        df.sort_values(by=['smiles_id', 'mean_score'], inplace=True, ascending=[True, False])

        # Save the ranking to a CSV file
        df.to_csv(f'{dest}_{top_k}.csv', index=False, float_format='%.4f')

    print(f'\nResults:', file=log)
    print(f'\tIdentified {len(reference_dataset.uniprot.unique())} potential targets based on the mean similarity scores of the top k reference compounds for {len(reference_dataset.smiles_id.unique())} SMILES strings.', file=log)
    print(f"\tAdditional details are stored in {dest}_{{k}}.csv", file=log)

    print(f'Retrieving the potential targets took {time.time() - time_s:.4f} seconds.', file=log)
