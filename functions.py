#############
# FUNCTIONS #
#############


# Import libraries ################################################
import os
# Current working directory
curr_dir = os.getcwd()
data_dir = os.path.join('data')
notebook_dir = os.path.join('notebooks')
project_dir = os.path.join('..')

import pandas as pd
import numpy as np
from chembl_webresource_client.new_client import new_client


##################################################################
gene_name = input('Enter gene name: ')

def get_bioactivity_data(gene_name):
    target = new_client.target

    # Searching for target in ChEMBL database
    target_query = target.search(gene_name)
    targets = pd.DataFrame.from_dict(target_query) # list of all targets in ChEMBL database

    # Select and retrieve bioactivity data caspase-2 (1st entry)
    selected_target = targets.target_chembl_id[0]

    # Retrieve bioactivity data for caspase-2
    activity = new_client.activity
    res = activity.filter(target_chembl_id=selected_target).filter(standard_type="IC50")
    # Converting to dataframe
    res_df = pd.DataFrame.from_dict(res)
    # Saving to csv
    if not os.path.exists(os.path.join(project_dir, data_dir)): os.mkdir(os.path.join(project_dir, data_dir)) # Creating a data directory, if it doesn't exist
    res_df.to_csv(os.path.join(project_dir, data_dir, f'{gene_name}_bioactivity_data_raw.csv'), index=False)

    # Dropping compounds without standard value
    df = res_df.dropna(subset=['standard_value'])
    
    # Lebeling compounds as either being active, intermediate or inactive
    #- active (IC50 < 1000 nm)
    #- intermediate (1000, 10000 nm)
    #- inactive (> 10000 nm)
    labeling = lambda x: 'active' if float(x) <= 1000 else 'intermediate' if float(x) <= 10000 else 'inactive'
    df['bioavtivity_class'] = df.standard_value.apply(labeling)
    # Keeping only the relevant columns
    cols_to_keep = ['molecule_chembl_id', 'canonical_smiles', 'standard_value', 'bioavtivity_class']
    df = df[cols_to_keep]
    # Saving final dataframe to csv
    df.to_csv(os.path.join(project_dir, data_dir, f'{gene_name}_bioactivity_data_preprocessed.csv'), index=False)
    
    return df


# def pre_process(df):
#     """
#     1. Labeling compounds as either being active, intermediate or inactive
#     - intermediate (1000, 10000 nm)
#     - inactive (> 10000 nm)
#     - inactive (> 10000 nm)
#     2. 
#     """

#     labeling = lambda x: 'active' if float(x) <= 1000 else 'intermediate' if float(x) <= 10000 else 'inactive'
#     df['bioavtivity_class'] = df.standard_value.apply(labeling)
#     # Keeping only the relevant columns
#     cols_to_keep = ['molecule_chembl_id', 'canonical_smiles', 'standard_value', 'bioavtivity_class']
#     df = df[cols_to_keep]
#     # Saving final dataframe to csv
#     df.to_csv(os.path.join(project_dir, data_dir, f'{gene_name}_bioactivity_data_preprocessed.csv'), index=False)

