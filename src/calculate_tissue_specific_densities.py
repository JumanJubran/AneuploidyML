import statistics as stats

import pandas as pd

from get_info import *

def get_ts_info(Dataset):
    '''
    This function returns two dictionaries, the first contains for each chromosome arm the number of tissue-specific genes for each tissue,
    the second contains for each chromosome arm the number of total genes.
    :param Dataset:
    :return: chromrosome_ts_info, chromosomal_total_genes
    '''
    chromosome_ts_info = {}
    chromosomal_total_genes = {}
    for rowNum, row in Dataset.iterrows():
        chromosome_arm = row['Chromosome_arm']
        row = row.to_dict()
        
        if chromosome_arm not in chromosome_ts_info:
            chromosome_ts_info[chromosome_arm] = {}
            chromosomal_total_genes[chromosome_arm] = 0
        chromosomal_total_genes[chromosome_arm] += 1
        
        for item in row:
            tissue, feature = split_name(item)
            if 'expression' == feature:
                
                if tissue not in chromosome_ts_info[chromosome_arm]:
                    chromosome_ts_info[chromosome_arm][tissue] = 0
                if row[item] == row[item] and row[item]!='':
                    if row[item] > 2:
                        chromosome_ts_info[chromosome_arm][tissue] += 1
    return chromosome_ts_info, chromosomal_total_genes

def median_per_tissue(chromrosome_ts_info):
    '''
    This function calculates the median number of tissue-specific genes since some tissues have sub-tissues in TRACE.
    :param chromrosome_ts_info:
    :return:
    '''
    chromosome_ts_info_median = {}

    for chromosome in chromrosome_ts_info:
        chromosome_ts_info_median[chromosome] = {}
        tissue_values = {}
        for tissue in chromrosome_ts_info[chromosome]:
            parent_tissue = tissue.split('-')[0].strip()    #Get the name of the parent_tissue, for example BRAIN that have sub-tissues now it will be called as BRAIn
            if parent_tissue not in tissue_values:
                tissue_values[parent_tissue] = []
            tissue_values[parent_tissue].append(chromrosome_ts_info[chromosome][tissue])

        for tissue in tissue_values:
            chromosome_ts_info_median[chromosome][tissue] = stats.median(tissue_values[tissue])
    return chromosome_ts_info_median



'''
    Function: calculate_ts_density
    Return: NONE!
    Input files: NONE!
    Input data:  chromrosome_ts_info_median, chromosomal_total_genes
    This function calculates the density of tissue-specific genes for each chromosome arm.
'''
def calculate_ts_density(chromrosome_ts_info_median, chromosomal_total_genes):
    ts_density_per_arm_tissue = {'Chromosome arm':[]}
    tissues = list(chromrosome_ts_info_median['1p'].keys())
    for t in tissues:
        ts_density_per_arm_tissue[t]=[]

    for chromosome in chromrosome_ts_info_median:
        ts_density_per_arm_tissue['Chromosome arm'].append(chromosome)

        for tissue in chromrosome_ts_info_median[chromosome]:
            ts_density_per_arm_tissue[tissue].append(chromrosome_ts_info_median[chromosome][tissue]/chromosomal_total_genes[chromosome])

    ts_density_per_arm_tissue_dataframe = pd.DataFrame.from_dict(ts_density_per_arm_tissue)
    ts_density_per_arm_tissue_dataframe.to_csv('Processed datasets/Tissue_specific_gene_densities.csv', index = False)
    return ts_density_per_arm_tissue_dataframe

def add_feature_to_dataset(Dataset, ts_gene_densities, data_type):
    '''
    This function adds the tissue specific densities to the model's features as a z-score per tissue.
    :param Dataset:
    :param ts_gene_densities:
    :return:
    '''
    ts_density_dict = {}
    cancer_tissue_dict = get_cancers_associated_tissue(data_type,'GTEx')                    #This function from get_info file
    density_feature = []
    for rowNum, row in ts_gene_densities.iterrows():
        row = row.to_dict()
        ts_density_dict[row['Chromosome arm']] = row
    
    for rowNum, row in Dataset.iterrows():
        disease_tissue  = []
        unaffected_tissue = []
        row = row.to_dict()
        cancer_type = row['Type']
        c_arm = row['Arm']
        
        d_tissue = cancer_tissue_dict[cancer_type]
        densities = ts_density_dict[c_arm]
        for tissue in densities:
            tissue_ = tissue.replace('_',' ').strip().lower()
            if tissue != 'Chromosome arm':
                if d_tissue in tissue_ and d_tissue!='':
                    disease_tissue.append(densities[tissue])
                else:
                    unaffected_tissue.append(densities[tissue])

        if len(disease_tissue) != 0:
            
            '''
            Z-score calculation
            '''
            chromosome_score = stats.median(disease_tissue)
            mean = stats.mean(disease_tissue+unaffected_tissue)
            stdv = stats.stdev(disease_tissue+unaffected_tissue)
            z_score = (chromosome_score - mean)/stdv
            density_feature.append(z_score)
            
        else:
            density_feature.append('')
        
    Dataset['tissue specific genes density'] = density_feature
    return Dataset


def ts_densities_main(Dataset_genes,arm_cancer_dataset, data_type):

    if data_type == 'TCGA':
        chromosome_ts_info, chromosomal_total_genes = get_ts_info(Dataset_genes)
        chromrosome_ts_info_median = median_per_tissue(chromosome_ts_info)
        ts_gene_densities = calculate_ts_density(chromrosome_ts_info_median, chromosomal_total_genes)
    else:
        ts_gene_densities = pd.read_csv('Processed datasets/Tissue_specific_gene_densities.csv')
    Dataset = add_feature_to_dataset(arm_cancer_dataset, ts_gene_densities, data_type)
    return Dataset
    
