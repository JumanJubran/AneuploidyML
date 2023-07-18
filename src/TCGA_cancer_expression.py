import statistics
import os
from get_info import *

def get_neutral_arm_per_sample(arm_level):
    chromosome_arms = get_chromosome_arms()     #from get_info
    # finding neutral chromosome arms per sample
    sample_neutral_arms = {}
    for rowNum, row in arm_level.iterrows():
        row = row.to_dict()
        sample_neutral_arms[row['Sample']] = []
        for item in row:
            if item in chromosome_arms:
                # If we want to calculate non copy number altered samples (neutral) we add those with zero value.
                if row[item] == 0:
                    sample_neutral_arms[row['Sample']].append(item)
    return sample_neutral_arms

def get_gene_neutral_samples_expression(sample_neutral_arms, row, chrom_arm):
    neutral_samples_expression = []
    for sample in row:

        if sample[:-1] in sample_neutral_arms:

            if chrom_arm in sample_neutral_arms[sample[:-1]]:
                neutral_samples_expression.append(row[sample])
    return neutral_samples_expression

def add_TCGA_features_to_dataset(Dataset, added_TCGA_features_to_dataset, gene_median_expression):
    for rowNum, row in Dataset.iterrows():
        gene_id = row['Gene_ID']
        if gene_id in gene_median_expression:
            for cancer in added_TCGA_features_to_dataset:
                added_TCGA_features_to_dataset[cancer].append(gene_median_expression[gene_id][cancer])
        else:
            for cancer in added_TCGA_features_to_dataset:
                added_TCGA_features_to_dataset[cancer].append('')
    for cancer in added_TCGA_features_to_dataset:
        Dataset[cancer.lower()+'_(neutarl) TCGA expression'] = added_TCGA_features_to_dataset[cancer]
    return Dataset

def gene_median_cancer_expression(Dataset):
    arm_level = pd.read_csv("../Datasets/Taylor_et_al._Arm-Level_WGD_TCGA_data.csv")

    sample_neutral_arms = get_neutral_arm_per_sample(arm_level)
    gene_chromosome, chromosome_genes = get_genes_id_chromosome_arm()      #from get_info

    gene_median_expression = {}
    flag_int = 0

    added_TCGA_features_to_dataset = {}
    files_array = os.listdir("../Datasets/TCGA expression/")
    
    for file in files_array:

        cancer_name = file.split('.')[0]
        
        added_TCGA_features_to_dataset[cancer_name] = []

        cancer_expression = pd.read_csv("../Datasets/TCGA expression/"+file, sep='\t')

        for rowNum, row in cancer_expression.iterrows():
            row = row.to_dict()

            if row['Ensembl_ID'][:15] in gene_chromosome:
                gene_id = row['Ensembl_ID'][:15]
                del row['Ensembl_ID']
                chrom = gene_chromosome[gene_id]
                neutral_samples_expression = get_gene_neutral_samples_expression(sample_neutral_arms, row, chrom)


                
                if flag_int == 0:
                    gene_median_expression[gene_id] = {}
                    if len(neutral_samples_expression)!=0:
                        gene_median_expression[gene_id][cancer_name] = statistics.median(neutral_samples_expression)
                    else:
                        gene_median_expression[gene_id][cancer_name] = ''
                
                else:
                    if len(neutral_samples_expression)!=0:
                        gene_median_expression[gene_id][cancer_name] = statistics.median(neutral_samples_expression)
                    else:
                        gene_median_expression[gene_id][cancer_name] = ''
        flag_int = 1

    Dataset = add_TCGA_features_to_dataset(Dataset, added_TCGA_features_to_dataset, gene_median_expression)
    return Dataset

def TCGA_Main(Trace_Dataset):
    Dataset = gene_median_cancer_expression(Trace_Dataset)
    Dataset.to_csv('Processed datasets/3_TRACE_Dataset_with_TCGA.csv', index=False)
    return Dataset
