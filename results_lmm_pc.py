## results model meixi 
## IMPORT ALLELES FREQ 
import numpy as np
import pandas as pd
import os
import glob
import matplotlib.pyplot as plt
import seaborn as sns

dict_offset_nooffset_partitions = pd.read_csv('/home/tbellagio/scratch/slim_grenenet/data/dict_offset_nooffset_partitions.csv')

pattern = os.path.join('/home/tbellagio/scratch/slim_grenenet/results', '**', 'lmm_pc_results10env.csv')
 
lmm_pc_results_files = glob.glob(pattern, recursive=True)

def add_missing_rows(results_lmm):
    max_value = results_lmm.index.str.split('.').str[1].astype(int).max()
    expected_rows = [f'result.{i}' for i in range(1,max_value + 1)]
    missing_rows = set(expected_rows) - set(results_lmm.index)
    #check 
    #set(results_lmm.index) -  set(expected_rows)

    # Initialize an empty list to store dictionaries
    list_of_dicts = []

    # Create dictionaries for each index value and add them to the list
    for index in missing_rows:
        new_row = {'index': index, 'R2m': np.nan, 'R2c': np.nan, 'beta': np.nan, 'beta_p': np.nan, 'BIC': np.nan}
        list_of_dicts.append(new_row)

    # Create a DataFrame from the list of dictionaries
    missing_rows = pd.DataFrame(list_of_dicts)

    results_lmm = results_lmm.reset_index()
    results_lmm = pd.concat([results_lmm,missing_rows])
    results_lmm['index'] = results_lmm['index'].str.split('.').str[1].astype(int)

    results_lmm = results_lmm.sort_values(by = 'index')

    results_lmm = results_lmm.reset_index(drop=True).drop('index',axis=1)
    #len(results_lmm)

    return results_lmm

def metrics_wtreshold(name, results_lmm, th, min_pvalue):

    results_lmm['sig'] = results_lmm['p_value_env']<= th

    true_positives = len(results_lmm[(results_lmm['is_within_range'] ==True) & (results_lmm['sig'] ==True)])

    true_negatives = len(results_lmm[(results_lmm['is_within_range'] ==False) & (results_lmm['sig'] ==False)])

    false_negatives = len(results_lmm[(results_lmm['is_within_range'] ==True) & (results_lmm['sig'] ==False)])

    false_positives = len(results_lmm[(results_lmm['is_within_range'] ==False) & (results_lmm['sig'] ==True)])

    if false_positives + true_positives == 0:
        fdr = np.nan  # Set FDR to 0 if the denominator is zero
    else:
        fdr = false_positives / (false_positives + true_positives)

    new_row = {'name': name,
               'true_positives': true_positives,
               'true_negatives': true_negatives,
               'false_negatives': false_negatives,
               'false_positives': false_positives,
               'fdr': fdr,
               'min_p_value': min_pvalue, 
               'th': str(th)}

    return new_row

result_metrics = pd.DataFrame(columns=['name', 'true_positives', 'true_negatives', 'false_negatives', 'false_positives', 'fdr', 'min_p_value', 'th'])

for result in lmm_pc_results_files: 

    name = result.split('/')[-3] + '_' + result.split('/')[-4] + '_' + result.split('/')[-5]

    allele_freq_norm_file = result.replace('/lmm/lmm_pc_results10env.csv','/allele_freq_norm10env.csv')
    causal_loci_file = result.split('arq')[0] + 'arq' + result.split('arq')[1].split('/')[0] + '/loci_effectsize.csv'

    results_lmm = pd.read_csv(result,index_col=[0])

    bc = 0.05 / len(results_lmm)

    ## first i need to check that all results are there 
    results_lmm = add_missing_rows(results_lmm)

    ## read allele freq norm to get chrom pos 
    allele_freq_norm = pd.read_csv(allele_freq_norm_file)
    results_lmm = pd.concat([allele_freq_norm['chrom_pos'], results_lmm],axis=1)

    #get the causal loci 
    causal_loci = pd.read_csv(causal_loci_file)
    causal_loci = causal_loci.merge(dict_offset_nooffset_partitions, left_on = 'pos', right_on = 'offset')
    results_lmm = results_lmm.merge(dict_offset_nooffset_partitions[['offset', 'partition']], left_on = 'chrom_pos', right_on = 'offset')
    results_lmm.loc[:, 'is_within_range'] = results_lmm['partition'].isin(causal_loci['partition'])


    ## get the min p value for teh records 
    min_pvalue = results_lmm['p_value_env'].min()


    for i in [0.05, 0.0005, 0.000005, 0.00000005, bc]:
        new_row = metrics_wtreshold(name, results_lmm, i, min_pvalue)
        result_metrics.loc[len(result_metrics)] = new_row
result_metrics.to_csv('results_lmm_pc_10env.csv')