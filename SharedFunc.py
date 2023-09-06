"""

This file is part of the Compound-Target Activity (CTA) Prediction Program.

    
Prepared By: Alhasbary
Date: 05/09/2023

"""
import argparse
import gc

import numpy as np
import pandas as pd
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

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


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
    print(f"    Number of unique chembl targets  = {len(df.target_chembl_id.unique())}", file=log)
    print(f"    Number of unique uniprot targets = {len(df.uniprot.unique())}", file=log)
    print(f"    Number of unique compounds       = {len(df.molecule_chembl_id.unique())}", file=log)
    print(f"    Target types included: {df.target_type.unique()}", file=log)
    df_a = df.groupby(['uniprot']).count()
    print('\nCompound statistics in the targets:', file=log)
    print('min = ', df_a.molecule_chembl_id.min(), '\t median = ', df_a.molecule_chembl_id.median(), '\t mean = ', df_a.molecule_chembl_id.mean(),
          '\t max = ', df_a.molecule_chembl_id.max(), file=log)


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
    desc = f'''1- Extracted bioactivity data from ChEMBL with the following selection constraints:
    - Canonical smiles (not null)
    - Accession (not null)
    - Assay type (binding assays)
    - Confidence score >= {args.confidence_score}
    - Standard units (nM)
    - Standard value <= {args.standard_value}
    - Standard relation (=)
    - Standard type (IC50, Ki, Kd, or EC50)
    - pChembl value (not null)
    - Data validity comment (null or manually validated)
    - Activity comment (not undetermined, inconclusive, unspecified, inactive, or not active)
    - Potential duplicate (0)
    
    Run the code with -h args to choose the desired settings.
    '''
    print(desc, file=log)

    sqlstr = '''select DISTINCT
        md.chembl_id as molecule_chembl_id, cs.canonical_smiles, act.standard_value, 
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
        ass.confidence_score in '''
    
    if type(args.confidence_score) == list:
        if len(args.confidence_score) > 1:
            sqlstr += str(tuple(args.confidence_score))
        else:
            sqlstr += '(' + str(args.confidence_score[0]) + ')'
    else:
        sqlstr += '(' + str(args.confidence_score) + ')'
    
    sqlstr += f" and act.standard_value <= {args.standard_value}"

    dbname = os.path.join(args.inputPath, "mini_chembl.db")
    try:
        time_s = time.time()
        selected_bioactivity = from_db_to_pd(sqlstr, dbname)
        selected_bioactivity = selected_bioactivity.drop_duplicates(
            subset=selected_bioactivity.columns.difference(['protein_class_desc']))
        print_statistics(log, selected_bioactivity)
        elapsed_time = time.time() - time_s
        print(f'Extraction completed in {elapsed_time:.2f} seconds.', file=log)
        return selected_bioactivity

    except NameError as e:
        print(f"An error occurred: {e}")
        exit(0)


def processMS(smi, args):
    """
    Process Molecular Structures.

    :param smi: SMILES string.
    :param args: argparse object containing program arguments.
    :return: List of SMILES strings with standardized structures. Problematic SMILES strings are replaced with 'nan' and removed.
    """  
    
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
    else:
        smiles = ''
    if smiles == '':
        smiles = np.nan
    return smiles

    
def preprocess_chembl(log, args, df_db):
    """
    Preprocess ChEMBL Molecular Structures.

    :param log: Log file for recording processing information.
    :param args: argparse object containing program arguments.
    :param df_db: DataFrame containing the selected ChEMBL data.
    :return: DataFrame containing the cleaned selected ChEMBL data.
    """

    desc = f'''2- Clean SMILES strings, preprocess them using RDKit and MolVS. The preprocessing steps include normalizing
            structures, desalting, and neutralizing charge. Users can choose to disable desalting and charge neutralization
            if desired. Additionally, the user can specify the number of CPU cores used for processing.
    
    Run the code with -h args to choose the desired settings.
    '''
    print(desc, file=log)
    
    # Select smiles strings to preprocess.
    smiles = df_db.drop_duplicates(subset='canonical_smiles', keep="first")[
        ['molecule_chembl_id', 'canonical_smiles']]

    print('ChEMBL data before preprocessing:', file=log)
    print_statistics(log, df_db)

    # Start parallel preprocessing
    time_s = time.time()
    smi = smiles.canonical_smiles.unique()
    cleaned_smi = Parallel(n_jobs=int(args.n_jobs), backend="multiprocessing", batch_size=len(smi) // 4)(
        delayed(processMS)(smiles, args) for smiles in smi)

    # Replace the old smiles with the preprocessed smiles.
    smiles.loc[:, 'canonical_smiles'] = cleaned_smi
    df_db.drop(labels='canonical_smiles', axis="columns", inplace=True)
    #df_db.reset_index(level=0, inplace=True, drop=True)
    df_db = pd.merge(smiles, df_db, on='molecule_chembl_id')

    # Remove rows containing problematic smiles
    cleaned_data = df_db[df_db['canonical_smiles'].notna()]
    cleaned_data.rename(columns={'canonical_smiles': 'smiles'})
    print('ChEMBL data after preprocessing:', file=log)

    # Print statistics and save the results
    print_statistics(log, cleaned_data)
    elapsed_time = time.time() - time_s
    print(f'ChEMBL preprocessing completed in {elapsed_time:.2f} seconds.', file=log)
    
    return cleaned_data


def preprocess_smiles(log, args, df_db):
    """
    Process NPs Molecular Structures.

    :param log: Log file for recording processing information.
    :param args: argparse object containing program arguments.
    :param df_db: DataFrame containing the NPs data.
    :return: DataFrame containing the cleaned NPs data.
    """

    # Select smiles strings to preprocess.
    df_db = df_db.drop_duplicates(subset='smiles', keep="first")
    #df_db.reset_index(level=0, inplace=True, drop=True)

    print('SMILES strings before preprocessing:', file=log)
    print(f'Number of compounds: {len(df_db.smiles)}', file=log)
    print(f'Number of unique compounds: {len(df_db.smiles.unique())}', file=log)

    # Start clean preprocess in parallel
    time_s = time.time()
    smi = df_db.smiles.unique()
    cleaned_smi = Parallel(n_jobs=int(args.n_jobs), backend="multiprocessing", batch_size=len(smi) // 4)(
        delayed(processMS)(smiles, args) for smiles in smi)

    # Replace the old smiles with the preprocessed smiles.
    df_db.drop(labels='smiles', axis="columns", inplace=True)
    df_db = df_db.copy()  # Create a copy of the DataFrame
    df_db['smiles'] = cleaned_smi

    # Remove any rows which contained problematic smiles
    cleaned_data = df_db[df_db['smiles'].notna()]


    # Print statistics and save the results
    print('SMILES strings after preprocessing:', file=log)
    print(f'Number of compounds: {len(cleaned_data.smiles)}', file=log)
    print(f'Number of unique compounds: {len(cleaned_data.smiles.unique())}', file=log)
    elapsed_time = time.time() - time_s
    print(f' preprocessing completed in {elapsed_time:.2f} seconds.', file=log)
    
    return cleaned_data


def chembl_fingerprints(log, args, cleaned_data):
    """
    Generating fingerprints for ChEMBL Molecular Structures.

    :param log: Log file for recording processing information.
    :param args: argparse object containing program arguments.
    :param cleaned_data: DataFrame containing the cleaned ChEMBL data.
    :return: list of valid fingerprints and the crosponding compound ids.
    """
    
    smiles = cleaned_data.drop_duplicates(subset='canonical_smiles', keep="first")[
        ['molecule_chembl_id', 'canonical_smiles']]

    print(f'\nStart generating fingerprints for {len(smiles)} ChEMBL compounds', file=log)
    smi = smiles.canonical_smiles
    df_compound_ids = list(smiles.molecule_chembl_id.values)
    valid_fps, compound_ids = gen_fps(log, args, smi, df_compound_ids)
    return valid_fps, compound_ids


def nps_fingerprints(log, args, cleaned_data):
    """
    Generating fingerprints for NPs Molecular Structures.

    :param log: Log file for recording processing information.
    :param args: argparse object containing program arguments.
    :param cleaned_data: DataFrame containing the cleaned NPs data.
    :return: list of valid fingerprints and the crosponding compound ids.
    """
    cleaned_data = cleaned_data.drop_duplicates(subset='smiles', keep="first")
    print(f'\nStart generating fingerprints for {len(cleaned_data)} smiles', file=log)
    smi = cleaned_data.smiles
    df_compound_ids = list(cleaned_data.smiles_id.values)
    valid_fps, compound_ids = gen_fps(log, args, smi, df_compound_ids)
    return valid_fps, compound_ids


def gen_fp(smi, args):
    mol = MolFromSmiles(smi)
    fp = 'Issue'
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
    return fp


def gen_fps(log, args, smi, compound_ids):
    time_s = time.time()
    fps = Parallel(n_jobs=int(args.n_jobs), backend="multiprocessing")(
        delayed(gen_fp)(smiles, args) for smiles in smi)

    invalid_fps_ids = list()
    valid_fps = list()
    j = 0
    for i in range(len(fps)):
        if fps[i] == 'Issue':
            invalid_fps_ids.append(compound_ids.pop(i - j))
            j += 1
        else:
            valid_fps.append(fps[i])
    print("invalid fingerprints ", len(invalid_fps_ids), file=log)
    print(f'Generating fingerprints for {len(valid_fps)} smiles took {time.time() - time_s: .4f} seconds.',
          file=log)
    return valid_fps, compound_ids


def ss(i, fp, fps, np_ids, chembl_ids, Tc=0.85):
    Tcs = list()
    ss_results = np.array(DataStructs.BulkTanimotoSimilarity(fp, fps))
    idx = np.where(ss_results >= Tc)[0]
    for j in idx:
        Tcs.append([np_ids[i], chembl_ids[j], ss_results[j]])
    if len(Tcs) > 0:
        del ss_results, idx
        gc.collect()
        return Tcs


def print_ss_results(log, args, results):
    df = results.copy()
    print(f"Results = {df.shape}", file=log)
    print("\n=====================================", file=log)
    print("Summary\t\tChEMBL\t\tNPs", file=log)
    print("=====================================", file=log)
    print(f"100 identical    = {len(df[df['ss'] == 1.0].molecule_chembl_id.unique())}",
          end='\t', file=log)
    print(f"{len(df[df['ss'] == 1.0].smiles_id.unique())}", file=log)

    print(f"[90, 99] similar = {len(df[(df['ss'] < 1.0) & (df['ss'] >= 0.9)].molecule_chembl_id.unique())}",
          end='\t', file=log)
    print(f"{len(df[(df['ss'] < 1.0) & (df['ss'] >= 0.9)].smiles_id.unique())}", file=log)

    print(f"[85, 90] similar = {len(df[(df['ss'] < 0.9) & (df['ss'] >= 0.85)].molecule_chembl_id.unique())}",
          end='\t', file=log)
    print(f"{len(df[(df['ss'] < 0.9) & (df['ss'] >= 0.85)].smiles_id.unique())}", file=log)

    print(f"Total Uniques    = {len(df.molecule_chembl_id.unique())} \t {len(df.smiles_id.unique())}", file=log)


def conductSS(log, args, reference_fps, reference_ids, smiles_fps, smiles_ids):
    """
    Conduct Chemical Similarity Searching.

    :param log: Log file for recording search information.
    :param args: argparse object containing program arguments.
    :param reference_fps: List of fingerprints for reference compound structures.
    :param reference_ids: List of IDs for reference compound structures.
    :param smiles_fps: List of fingerprints for query compound structures.
    :param smiles_ids: List of IDs for query compound structures.
    :return: DataFrame containing the results of chemical similarity searching.
    """
    desc = f'''3- Perform chemical similarity searches between two datasets. RDKit is used to generate molecular fingerprints 
    for each compound. The cleaned SMILES strings for each compound are converted into {args.nBits}-bit {args.fingerprint} fingerprints 
    with a radius of {args.radius}. Targets that interact with at least one compound having a Tanimoto coefficient (Tc) of 
    {args.Tc_score} or higher with any query compound are considered and saved.

    To generate fingerprints for SMILES strings, the user can choose one of four specified fingerprint types: Avalon, 
    Morgan ECFP, Morgan FCFP, and MACCS. The user can also specify the desired settings for nBits to control the length 
    of the generated fingerprint and the radius parameter, which governs the size of the substructures used to generate 
    the fingerprint. Additionally, the user can set a custom Tc similarity threshold.

    The code utilizes the Parallel function from the joblib library to perform similarity searches with a batch size of 
    {args.batch}. This allows for efficient processing while addressing memory constraints. Users can specify the desired 
    batch size.

    Run the code with -h args to choose the desired settings.
    '''
    print(desc, file=log)

    # Start preprocessing in chunks for memory issues
    print(f'\nStart chemical similarity searches for {len(smiles_fps)} smiles\n', file=log)
    time_s = time.time()
    result = list()
    if len(smiles_fps) >= 1000:
        results = Parallel(n_jobs=int(args.n_jobs), backend="multiprocessing", batch_size=args.batch)(
            delayed(ss)(i, fp, reference_fps, smiles_ids, reference_ids, args.Tc_score) for i, fp in enumerate(smiles_fps))

        for x in results:
            if x is not None:
                for s in x:
                    result.append(s)
    else:
        Tc_scores = np.zeros(shape=(len(smiles_fps), len(reference_ids)), dtype=np.float16)
        for i, fp in enumerate(smiles_fps):
            Tc_scores[i] = DataStructs.BulkTanimotoSimilarity(fp, reference_fps)
        for i in range(Tc_scores.shape[0]):
            for j in range(Tc_scores.shape[1]):
                if Tc_scores[i][j] >= args.Tc_score:
                    result.append([smiles_ids[i], reference_ids[j], Tc_scores[i][j]])

    # Save the results

    if len(result) > 0:
        df_result = pd.DataFrame(result, columns=['smiles_id', 'molecule_chembl_id', 'ss'])
        df_result['ss'] = df_result['ss'].astype(float)        
        print_ss_results(log, args, df_result)
        elapsed_time = time.time() - time_s
        print(f' Simlarity searching completed in {elapsed_time:.2f} seconds.', file=log)        
        return df_result
    else:
        print('No similarities were found!', file=log)
    

def potential_targets_ss(log, args, cleaned_ChEMBL, ChEMBL_similarity, dest):
    """
    Save the results generated by conductSS.

    :param log: Log file for recording search information.
    :param args: argparse object containing program arguments.
    :param cleaned_ChEMBL: DataFrame containing the ChEMBL dataset with target information.
    :param ChEMBL_similarity: DataFrame containing the results of conductSS.
    :param dest: File name for saving the potential targets for each SMILES string.
    """
    desc = f'''4- Retrieve potential targets for each SMILES string based on similarity searching.\n'''
    print(desc, file=log)
    
    selected_ChEMBL = cleaned_ChEMBL.copy()
    similarities = ChEMBL_similarity.copy()

    # Consider a similarity threshold of {args.Tc_score} or better.
    chemble_id_with_Tc_similarity = similarities[similarities['ss'] >= args.Tc_score].molecule_chembl_id.unique()
    PT_with_Tc_similarity = selected_ChEMBL[selected_ChEMBL.molecule_chembl_id.isin(chemble_id_with_Tc_similarity)]

    # Generate a compound-target activity pair dataset with a similarity threshold of {args.Tc_score} or better.
    PT_ss = pd.merge(similarities, PT_with_Tc_similarity, on='molecule_chembl_id')
    print(f'Chemical similarity searches with {args.Tc_score} similarity threshold results:', file=log)
    print(f'\t{len(PT_ss.uniprot.unique())} potential targets for {len(PT_ss.smiles_id.unique())} SMILES strings.', file=log)
    print(f"More details are stored in {dest}", file=log)
    
    # Save related information
    columns = ['smiles_id', 'molecule_chembl_id', 'ss', 'standard_type', 'pchembl_value', 'target_chembl_id', 'uniprot']
    PT_ss = PT_ss[columns]
    PT_ss = PT_ss.drop_duplicates()
    #PT_ss.reset_index(level=0, inplace=True, drop=True)
    PT_ss.sort_values(by=['smiles_id', 'ss'], inplace=True, ascending=[True, False])
    PT_ss.to_csv(dest, index=False, float_format='%.4f')


def createCTA(log, args, cleaned_ChEMBL, ChEMBL_similarity, dest):
    """
    Creating a Compound-Target Activity (CTA) Pair Dataset.

    :param log: Log file for recording dataset creation information.
    :param args: argparse object containing program arguments.
    :param cleaned_ChEMBL: DataFrame containing the ChEMBL dataset with target information.
    :param ChEMBL_similarity: DataFrame containing the results of chemical similarity searching.
    :param dest: File name for saving the CTA dataset.
    """
    
    desc = f'''5- Create a compound-target pair dataset. Only compounds with a similarity score of {args.Tc_score} or 
    higher are considered from the cleaned ChEMBL dataset. When identical target-compound pairs are found, the {args.agg} 
    activity value is used. Additionally, the dataset is pruned, retaining only compound-target pairs where targets have 
    at least {args.minCompounds} compounds. Furthermore, the number of compounds per target is capped at {args.maxCompounds} 
    through random sampling. 
    Run the code with -h args to choose the desired settings.\n'''
    print(desc)
    time_s = time.time()
    selected_ChEMBL = cleaned_ChEMBL.copy()
    similarities = ChEMBL_similarity.copy()

    # Consider compounds with a similarity score of {args.Tc_score} or higher.
    chemble_id_with_Tc_similarity = similarities[similarities['ss'] >= args.Tc_score].molecule_chembl_id.unique()
    CTA_with_Tc_similarity = selected_ChEMBL[selected_ChEMBL.molecule_chembl_id.isin(chemble_id_with_Tc_similarity)]

    # Create a compound-target activity pair dataset considering {args.Tc_score} or better similarity.
    CTA_pair = selected_ChEMBL[selected_ChEMBL.target_chembl_id.isin(CTA_with_Tc_similarity.target_chembl_id)]

    # When identical target-compound pairs are found, the {args.agg} activity value is used.
    column_map = {col: "first" for col in CTA_pair.columns}
    column_map["pchembl_value"] = args.agg
    CTA_multiple_activity = CTA_pair.groupby(["uniprot", "standard_type", "canonical_smiles"],
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
    print('After preprocessing the data for machine learning models,', file=log)
    print_statistics(log, CTA)
    
    # Save related information
    columns = ['molecule_chembl_id', 'standard_type', 'pchembl_value', 'target_chembl_id', 'uniprot']
    CTA = CTA[columns]
    CTA = CTA.drop_duplicates()
    #CTA.reset_index(level=0, inplace=True, drop=True)
    CTA.sort_values(by=['molecule_chembl_id', 'pchembl_value'], inplace=True, ascending=[True, False])    
    CTA.to_csv(dest, index=False, float_format='%.4f')
    
    print(f'\nCreating the compound-target activity pair dataset took {time.time() - time_s:.4f} seconds.', file=log)

