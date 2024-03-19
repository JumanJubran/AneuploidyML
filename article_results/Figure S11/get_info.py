

import pandas as pd


'''
    Function: get_cancers_associated_tissue
    Return: cancer_tissue_dict
    Input files: cancer tissue association.csv
    Input data: NONE!
    This function returns a dictionary of cancer types and their disease tissue.
'''
def get_cancers_associated_tissue(data_type,feature_type):
    
    if data_type == 'TCGA':
        cancer_tissue = pd.read_csv("../../Datasets/cancer tissue association.csv")
    if data_type == 'CCL':
        cancer_tissue = pd.read_csv("../Datasets/CCL data/cancer_tissue_association.csv")
    cancer_tissue_dict = {}
    
    if feature_type == 'DepMap':
        tissue = 'Tissue - DepMap'
    else:
        tissue = 'Tissue - GTEx'
    for rowNum, row in cancer_tissue.iterrows():
        if row[tissue] == row[tissue]:
            cancer_tissue_dict[row['Type']] = row[tissue].lower()
    
    return cancer_tissue_dict
    

'''
    Function: get_genes_chromosome_arm
    Return: gene_chromosome, chromosome_genes
    Input files: file contain for each gene the chromosome arm its located on.
    Input data: NONE!
    This function returns two dictionaries, the first contains for each gene its location on a specific chromosome arm, 
    the second contains for each chromosome arm the genes that are located on it.
'''
def get_genes_chromosome_arm():
    chromosome_genes = {}
    gene_chromosome = {}
    data = pd.read_csv('../../Datasets/protein coding genes - chromosomal location(bioMart).csv')
    for rowNum, row in data.iterrows():
        chromosome = str(row['Chromosome/scaffold name']) + str(row['Karyotype band'])[0]
        gene_chromosome[row['Gene name']] = chromosome
        if chromosome not in chromosome_genes:
            chromosome_genes[chromosome] = []
        chromosome_genes[chromosome].append(row['Gene name'])
    return gene_chromosome, chromosome_genes

def get_genes_id_chromosome_arm():
    chromosome_genes = {}
    gene_chromosome = {}
    data = pd.read_csv("Processed datasets/1_TRACE_Dataset_processed_features.csv")
    for rowNum, row in data.iterrows():
        row = row.to_dict()
        gene_chromosome[row['Gene_ID']] = row['Chromosome_arm']
        if row['Chromosome_arm'] not in chromosome_genes:
            chromosome_genes[row['Chromosome_arm']] = []
        chromosome_genes[row['Chromosome_arm']].append(row['Gene_ID'])
    return gene_chromosome, chromosome_genes



def get_genes_name():
    gene_id_to_name = {}
    gene_name_to_id = {}
    data = pd.read_csv("Processed datasets/1_TRACE_Dataset_processed_features.csv")
    for rowNum, row in data.iterrows():
        row = row.to_dict()
        gene_id_to_name[row['Gene_ID']] = row['Gene_name']
        gene_name_to_id[row['Gene_name']] = row['Gene_ID']
    return gene_id_to_name, gene_name_to_id
'''
    Function: add_labels
    Return: write_row
    Input files: NONE!
    Input data: fieldnames
    This function returns the dictionary of output file's fieldnames
'''
def add_labels(fieldnames):
    write_row = {}
    for f in fieldnames:
        write_row[f] = ''
    return write_row
    
 
 
'''
    Function: get_chromosome_arms
    Return: chromosome_arms
    Input files: NONE!
    Input data: NONE!
    This function returns the 39 chromosome arms that have genes.
'''
def get_chromosome_arms():
    chromosome_arms = []
    for i in range(1,23):
        chromosome_arms.append(str(i)+'p')
        chromosome_arms.append(str(i)+'q')

    chromosome_arms.remove('13p')
    chromosome_arms.remove('14p')
    chromosome_arms.remove('15p')    
    chromosome_arms.remove('21p')
    chromosome_arms.remove('22p')    
    return chromosome_arms
    
    

'''
    Function: get_cancer_types
    Return: NONE!
    Input files: NONE!
    Input data: NONE!
    This function returns the 33 cancer types.
'''
def get_cancer_types():
    return ['LUSC', 'PCPG', 'THCA', 'DLBC', 'THYM', 'UCS', 'KIRC', 'MESO', 'PAAD', 'CESC',
    'LGG', 'TGCT', 'OV', 'ACC', 'STAD', 'KIRP', 'KICH', 'HNSC', 'READ', 'SARC', 'LIHC', 'SKCM',
    'UVM', 'UCEC', 'ESCA', 'LUAD', 'LAML', 'BRCA', 'COAD', 'GBM', 'PRAD', 'BLCA', 'CHOL']
    
    
'''
    Function: get_cancer_types
    Return: NONE!
    Input files: NONE!
    Input data: data_type - either TCGA or CCL
    This function returns the 24 cancer types and 10 CCLs.
'''
def get_cancer_types_in_model(data_type):
    if data_type == 'TCGA':
        return ['LUSC', 'THYM', 'UCS', 'KIRC', 'PAAD', 'CESC',
        'LGG', 'TGCT', 'OV', 'ACC', 'STAD', 'KIRP', 'KICH',  'LIHC', 'SKCM',
        'UCEC', 'ESCA', 'LUAD', 'LAML', 'BRCA', 'COAD', 'GBM', 'PRAD', 'BLCA']
    else:
        return ['Haematopoietic_and_lymphoid', 'Lung', 'Ovary', 'Skin', 'Breast', 'Pancreas', 'Central_nervous_system', 'Stomach', 'Kidney', 'Liver']




def split_name(column_name):

    split_name = column_name.split('_')
    tissue = split_name[0]
    feature = ''
    for i in range(1,len(split_name)):
        feature += split_name[i]
        feature += '_'
    
    return tissue, feature[:-1]
    
    

def get_num_genes_chromosome():
    data = pd.read_csv("/gpfs0/estiyl/users/juman/Projects/Aneuploidy most important/Data/Gene_Base Dataset.csv")
    chrom_gene_num = {}
    for rowNum, row in data.iterrows():
        if row['Chromosome arm'] not in chrom_gene_num:
            chrom_gene_num[row['Chromosome arm']] = 0
        chrom_gene_num[row['Chromosome arm']] += 1
    return chrom_gene_num
