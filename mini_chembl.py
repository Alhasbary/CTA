
"""
Create a custom mini-ChEMBL SQLite database tailored to fulfill specific application requirements. 
In our case, it significantly reduces the storage size from 22.4 GB (ChEMBL32) to just 714 MB.

Part of the Compound-Target Activity (CTA) Prediction Program.

Inputs:
   - data: Full path to the ChEMBL dataset [Required].
   - destination: Full path to save the custom ChEMBL version as 'mini_chembl.db' [Optional].

Output: 
   - 'mini_chembl.db' is stored in the specified output directory, containing customized information for the application's use.

Usage:
   - For help:
     python mini_chembl.py -h
   - To run the program:
     python mini_chembl.py --data=FullpathToDatabaseFile --destination=FullPathToSavedDatabase

Usage Example:
    Download the ChEMBL SQLite dataset from https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/ 
    tar -xvzf chembl_32_sqlite.tar.gz 
    python mini_chembl.py --data=chembl_32/chembl_32_sqlite/chembl_32.db --destination=input
    
Prepared By: Alhasbary
Date: 05/09/2023
"""


import sqlite3
import time
import argparse
import os


def parse_args():
    description = """Create a custom mini-ChEMBL SQLite database tailored for specific application needs."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--data", action="store", dest='dataPath', type=str, default=None,
                        help='Full filepath of the SQLite file containing the ChEMBL database.')
    parser.add_argument("--destination", action="store", dest='destinationPath', type=str,
                        default='input', help='Full filepath for saving the custom ChEMBL database. Default is "input".')
    args = parser.parse_args()

    return args


def create_mini_chembl_database(dataPath, destinationPath='input'):
    """
    Create a custom mini-ChEMBL SQLite database tailored for specific application needs.

    :param 
        dataPath: Full filepath of the SQLite file containing the ChEMBL database.
        destinationPath: Full filepath for saving the custom ChEMBL version database. Default is 'input'.
    :return:
        A custom_mini_chembl_db: custom mini ChEMBL32 = 714 MB
        
    """
    try:
        # Connect to the ChEMBL database
        conn = sqlite3.connect(os.path.join("", dataPath))
        cursor = conn.cursor()
        
        # Perform custom operations to create the mini-ChEMBL database        
        mini_chembl = sqlite3.connect(os.path.join(destinationPath, "mini_chembl.db"))
        mini_chembl_cur = mini_chembl.cursor()

        # Make some fresh tables using executescript()
        mini_chembl_cur.executescript('''
        DROP TABLE IF EXISTS molecule_dictionary;
        DROP TABLE IF EXISTS compound_structures;
        DROP TABLE IF EXISTS compound_properties;
        DROP TABLE IF EXISTS component_sequences;
        DROP TABLE IF EXISTS protein_classification;
        DROP TABLE IF EXISTS component_class;
        DROP TABLE IF EXISTS activities;
        DROP TABLE IF EXISTS target_components;
        DROP TABLE IF EXISTS target_dictionary;
        DROP TABLE IF EXISTS assays;    

        CREATE TABLE molecule_dictionary (
        molregno BIGINT NOT NULL,
        chembl_id VARCHAR(20) NOT NULL,
        CONSTRAINT pk_moldict_molregno PRIMARY KEY (molregno),
        CONSTRAINT uk_moldict_chemblid UNIQUE (chembl_id)
        );

        CREATE TABLE compound_structures (
        molregno BIGINT NOT NULL,
        canonical_smiles VARCHAR(4000),
        CONSTRAINT pk_cmpdstr_molregno PRIMARY KEY (molregno),
        CONSTRAINT fk_cmpdstr_molregno FOREIGN KEY(molregno) REFERENCES molecule_dictionary (molregno) ON DELETE CASCADE
        );

        CREATE TABLE compound_properties (
        molregno BIGINT NOT NULL,
        np_likeness_score NUMERIC(3, 2),
        CONSTRAINT pk_cmpdprop_molregno PRIMARY KEY (molregno),
        CONSTRAINT fk_cmpdprop_molregno FOREIGN KEY(molregno) REFERENCES molecule_dictionary (molregno) ON DELETE CASCADE
        );

        CREATE TABLE component_sequences (
        component_id BIGINT NOT NULL,
        accession VARCHAR(25),
        sequence TEXT,     
        CONSTRAINT pk_targcomp_seqs_compid PRIMARY KEY (component_id),
        CONSTRAINT uk_targcomp_seqs_acc UNIQUE (accession)
        );
        
        CREATE TABLE protein_classification (
        protein_class_id BIGINT NOT NULL, 
        protein_class_desc VARCHAR(410) NOT NULL,
        CONSTRAINT prot_class_pk PRIMARY KEY (protein_class_id)
        ); 
        
        CREATE TABLE component_class (
        component_id BIGINT NOT NULL,
        protein_class_id BIGINT NOT NULL,
        comp_class_id BIGINT NOT NULL,
        CONSTRAINT pk_comp_class_id PRIMARY KEY (comp_class_id),
        CONSTRAINT fk_comp_class_compid FOREIGN KEY(component_id) REFERENCES component_sequences (component_id) ON DELETE CASCADE,
        CONSTRAINT fk_comp_class_pcid FOREIGN KEY(protein_class_id) REFERENCES protein_classification (protein_class_id) ON DELETE CASCADE,
        CONSTRAINT uk_comp_class UNIQUE (component_id, protein_class_id)
        );

        CREATE TABLE activities (
        activity_id BIGINT NOT NULL,
        assay_id BIGINT NOT NULL,
        molregno BIGINT,
        standard_relation VARCHAR(50),
        standard_value NUMERIC,
        standard_units VARCHAR(100),
        standard_type VARCHAR(250),
        activity_comment VARCHAR(4000),
        data_validity_comment VARCHAR(30),
        potential_duplicate SMALLINT,
        pchembl_value NUMERIC(4, 2),
        CONSTRAINT pk_act_activity_id PRIMARY KEY (activity_id),
        CONSTRAINT fk_act_assay_id FOREIGN KEY(assay_id) REFERENCES assays (assay_id) ON DELETE CASCADE,
        CONSTRAINT fk_act_molregno FOREIGN KEY(molregno) REFERENCES molecule_dictionary (molregno) ON DELETE CASCADE,
        CONSTRAINT ck_potential_dup CHECK (POTENTIAL_DUPLICATE IN (0,1)),
        CONSTRAINT ck_stand_relation CHECK (standard_relation in ('>','<','=','~','<=','>=','<<','>>'))
        );

        CREATE TABLE target_components (
        tid BIGINT NOT NULL,
        component_id BIGINT NOT NULL,
        targcomp_id BIGINT NOT NULL,
        CONSTRAINT pk_targcomp_id PRIMARY KEY (targcomp_id),
        CONSTRAINT fk_targcomp_compid FOREIGN KEY(component_id) REFERENCES component_sequences (component_id) ON DELETE CASCADE,
        CONSTRAINT fk_targcomp_tid FOREIGN KEY(tid) REFERENCES target_dictionary (tid) ON DELETE CASCADE,
        CONSTRAINT uk_targcomp_tid_compid UNIQUE (tid, component_id)
        );

        CREATE TABLE target_dictionary (
        tid BIGINT NOT NULL,
        target_type VARCHAR(30),
        pref_name VARCHAR(200) NOT NULL,
        organism VARCHAR(150),
        chembl_id VARCHAR(20) NOT NULL,
        CONSTRAINT pk_targdict_tid PRIMARY KEY (tid),
        CONSTRAINT fk_targdict_target_type FOREIGN KEY(target_type) REFERENCES target_type (target_type) ON DELETE CASCADE,
        CONSTRAINT uk_targdict_chemblid UNIQUE (chembl_id)
        );

        CREATE TABLE assays (
        assay_id BIGINT NOT NULL,
        description VARCHAR(4000),
        assay_type VARCHAR(1),
        tid BIGINT, 
        confidence_score SMALLINT,
        CONSTRAINT pk_assays_assay_id PRIMARY KEY (assay_id),
        CONSTRAINT fk_assays_tid FOREIGN KEY(tid) REFERENCES target_dictionary (tid) ON DELETE CASCADE
        );
        ''')

        # Copy molecule_dictionary table...
        sqlstr = 'SELECT molregno, chembl_id FROM molecule_dictionary'
        for row in cursor.execute(sqlstr):
            mini_chembl_cur.execute('''INSERT INTO molecule_dictionary (molregno, chembl_id)
            VALUES (?, ?)''', (row[0], row[1]))
        
        # Copy compound_structures table...
        sqlstr = 'SELECT molregno, canonical_smiles FROM compound_structures where canonical_smiles is not null'
        for row in cursor.execute(sqlstr):
            mini_chembl_cur.execute('''INSERT INTO compound_structures (molregno, canonical_smiles)
             VALUES (?, ?)''', (row[0], row[1]))
        
        # Copy compound_properties table...
        sqlstr = 'SELECT molregno, np_likeness_score FROM compound_properties'
        for row in cursor.execute(sqlstr):
            mini_chembl_cur.execute('''INSERT INTO compound_properties (molregno, np_likeness_score)
             VALUES (?, ?)''', (row[0], row[1]))
        
        # Copy component_sequences table...
        sqlstr = 'SELECT component_id, accession, sequence FROM component_sequences where accession is not null'
        for row in cursor.execute(sqlstr):
            mini_chembl_cur.execute('''INSERT INTO component_sequences (component_id, accession, sequence) 
            VALUES (?, ?, ?)''', (row[0], row[1], row[2]))
        
        # Copy protein_classification table...
        sqlstr = 'SELECT protein_class_id, protein_class_desc FROM protein_classification'
        for row in cursor.execute(sqlstr):
            mini_chembl_cur.execute('''INSERT INTO protein_classification (protein_class_id, protein_class_desc) 
            VALUES (?, ?)''', (row[0], row[1]))
        
        # Copy component_class table...
        sqlstr = 'SELECT component_id, protein_class_id, comp_class_id FROM component_class'
        for row in cursor.execute(sqlstr):
            mini_chembl_cur.execute('''INSERT INTO component_class (component_id, protein_class_id, comp_class_id) 
            VALUES (?, ?, ?)''', (row[0], row[1], row[2]))
        
        # Copy activities table...
        sqlstr = '''SELECT activity_id, assay_id, molregno, standard_relation, standard_type, standard_value, 
                 standard_units, activity_comment, data_validity_comment, potential_duplicate,pchembl_value 
                 FROM activities act WHERE
                            act.standard_value is not null and
                            act.pchembl_value is not null and
                            act.standard_units = 'nM' and
                            act.standard_relation = '=' and
                            act.potential_duplicate = 0 and
                            (
                            act.standard_type = 'IC50' or 
                            act.standard_type = 'Ki' or 
                            act.standard_type = 'Kd' or 
                            act.standard_type = 'EC50'
                            ) and
                            (
                            act.data_validity_comment is null or
                            act.data_validity_comment = 'manually validated'
                            )  and            
                            (
                            act.activity_comment is null or
                            lower(act.activity_comment) not in ('not active', 'inactive', 'inconclusive', 'undetermined') 
                            ) and
                            act.standard_value <= 10000'''
        for row in cursor.execute(sqlstr):
            mini_chembl_cur.execute('''INSERT INTO activities (activity_id, assay_id, molregno, standard_relation, 
            standard_type, standard_value,standard_units, activity_comment, data_validity_comment, potential_duplicate, 
            pchembl_value) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''', (row[0], row[1], row[2], row[3], row[4], row[5],
                                                                         row[6], row[7], row[8], row[9], row[10]))
                
        # Copy target_components table...
        sqlstr = 'SELECT targcomp_id, tid, component_id FROM target_components'
        for row in cursor.execute(sqlstr):
            mini_chembl_cur.execute('''INSERT INTO target_components (targcomp_id, tid, component_id) 
            VALUES (?, ?, ?)''', (row[0], row[1], row[2]))
        
        # Copy target_dictionary table...
        sqlstr = '''SELECT tid, target_type, pref_name, organism, chembl_id FROM target_dictionary
                    WHERE lower(target_type) IN ('single protein', 'protein complex')'''
        for row in cursor.execute(sqlstr):
            mini_chembl_cur.execute('''INSERT INTO target_dictionary (tid, target_type, pref_name, organism, chembl_id)
            VALUES (?, ?, ?, ?, ?)''', (row[0], row[1], row[2], row[3], row[4]))
        
        # Copy assays table...
        sqlstr = '''SELECT assay_id, assay_type, description, tid, confidence_score FROM assays
                    WHERE assay_type = 'B' '''
        for row in cursor.execute(sqlstr):
            mini_chembl_cur.execute('''INSERT INTO assays (assay_id, assay_type, description, tid, confidence_score) 
            VALUES (?, ?, ?, ?, ?)''', (row[0], row[1], row[2], row[3], row[4]))
        
       

        # Commit the changes
        mini_chembl.commit()
        conn.close()
        mini_chembl_cur.close()

        print("Custom mini-ChEMBL SQLite database created successfully!")

    except Exception as e:
        print(f"An error occurred: {e}")
        exit(0)


if __name__ == "__main__":
    args = parse_args()
    if args.dataPath:
        # Create a folder to save the custom ChEMBL database if it doesn't exist
        if not os.path.exists(args.destinationPath):
            os.makedirs(args.destinationPath)
        time_s = time.time()
        create_mini_chembl_database(args.dataPath, args.destinationPath)
        elapsed_time = time.time() - time_s
        print(f'Creation completed in {elapsed_time:.2f} seconds.')

    else:
        print("Please provide the path to the ChEMBL database using the --data argument.")
