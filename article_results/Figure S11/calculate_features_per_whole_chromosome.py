import pandas as pd
import statistics as stat
from calculate_tissue_specific_densities import *
from calculate_essential_densities import *
from essential_gene_density_Nichols import *
from standardize_expression_features import *


def fix_gene_zscore_file():
    data = pd.read_csv('../../src/Processed datasets/6_Genes z score.csv')
    chromosome = []
    for rowNum,row in data.iterrows():
        chromosome.append(row['Chromosome_arm'].replace('p','').replace('q',''))
    data.insert(loc=2, column='Chromosome',value = chromosome)
    data = data.drop(columns='Chromosome_arm')
    data.to_csv('processed_datasets/6_Genes z score.csv', index = False)


def calculate_median_top(Dataset):
    chromosome_hierarchical_info = {}  ## Building set of values based on chromosome and feature  dict[chromosome][feature] = [values]
    chromosome_hierarchical_gene_info = {}  ## Building set of values based on chromosome, feature and gene  dict[chromosome][feature][genes] = value
    chromosome_genes_count = {}  ## Number of genes on chromosome

    for rowNum, row in Dataset.iterrows():  ## Iteration on z-score dataset
        row = row.to_dict()

        ### Building the dataset above (chromosome_hierarchical_info,chromosome_hierarchical_gene_info,chromosome_genes_count)
        chromosome = row['Chromosome']
        if chromosome not in chromosome_hierarchical_info:
            chromosome_hierarchical_info[chromosome] = {}
            chromosome_hierarchical_gene_info[chromosome] = {}
            chromosome_genes_count[chromosome] = 0
            for feature in row:
                if feature not in ['Chromosome', 'Gene_name', 'Gene_ID']:
                    chromosome_hierarchical_info[chromosome][feature] = []
                    chromosome_hierarchical_gene_info[chromosome][feature] = {}

        chromosome_genes_count[chromosome] += 1
        for feature in row:
            if feature not in ['Chromosome', 'Gene_name', 'Gene_ID']:
                if row[feature] == row[feature] and row[feature] != '':
                    chromosome_hierarchical_info[chromosome][feature].append(row[feature])
                    chromosome_hierarchical_gene_info[chromosome][feature][row['Gene_name']] = row[feature]


    median_top_10_per = {}
    median_top_10_per['chromosome'] = []

    ## Calculating the median value of top 10% genes
    for chromosome in chromosome_hierarchical_info:
        if chromosome != '21p':
            median_top_10_per['chromosome'].append(chromosome)
            for feature in chromosome_hierarchical_info[chromosome]:

                if feature not in median_top_10_per:
                    median_top_10_per[feature] = []

                sorted_vals = sorted(chromosome_hierarchical_info[chromosome][feature])

                ten_precent = round(len(sorted_vals) * 0.1)

                if len(sorted_vals) != 0:
                    median_top_10_per[feature].append(stat.median(sorted_vals[-ten_precent:]))


                else:
                    median_top_10_per[feature].append('')

    pd.DataFrame.from_dict(median_top_10_per).to_csv('processed_datasets/7_median_top_10_percent.csv',index = False)

def get_features_types(Dataset):
    '''
    This functions extract the features (not including the tissue).
    :param Dataset:
    :return:
    '''
    features_types = []
    for item in Dataset.columns.tolist():
        if item not in ['chromosome arm']:
            item = item.split('_')
            curr_feature = ''
            for i in range(1,len(item)):
                curr_feature += item[i]
                curr_feature += ' '
            curr_feature = curr_feature[:-1]
            if curr_feature not in features_types:
                features_types.append(curr_feature)
    return features_types


def assign_fieldnames(features_types):
    fieldnames = {'Chromosome': [], 'Type': []}

    for f in features_types:
        fieldnames[f] = []

    return fieldnames


def get_features_values(row, g_tissue, e_tissue, c_type):
    features_values = {}
    for item in row:
        if item not in ['chromosome']:
            item2 = item.split('_')

            # Getting current tissue and curr feature

            curr_tissue = item2[0].strip().lower().replace(' ', '').replace('-', '')

            curr_feature = ''
            for i in range(1, len(item2)):
                curr_feature += item2[i]
                curr_feature += ' '
            curr_feature = curr_feature[:-1]

            if 'TCGA' in curr_feature:

                if c_type.replace(' ', '') in curr_tissue:

                    if curr_feature not in features_values:
                        features_values[curr_feature] = []

                    features_values[curr_feature].append(row[item])

            else:
                if 'rnai' in curr_feature or 'crispr' in curr_feature:
                    d_tissue = e_tissue
                else:
                    d_tissue = g_tissue

                if d_tissue in curr_tissue:
                    if curr_feature not in features_values:
                        features_values[curr_feature] = []

                    features_values[curr_feature].append(row[item])

    return features_values


def calculation_per_Chrom(z_score_dataset, disease_tissue_association, chrom_cancer_dataset):

    for rowNum,row in z_score_dataset.iterrows():

        row = row.to_dict()


        for rowNum2, row2 in disease_tissue_association.iterrows():

            chrom_cancer_dataset['Chromosome'].append(row['chromosome'])
            chrom_cancer_dataset['Type'].append(row2['Type'])

            g_tissue = row2['Tissue - GTEx'].replace(" ","").lower()


            if str(row2['Tissue - DepMap']) != 'nan':
                e_tissue = row2['Tissue - DepMap'].replace(" ","").lower()

            else:
                e_tissue = 'nan'

            features_values = get_features_values(row,g_tissue,e_tissue,row2['Type'].lower())



            for f in chrom_cancer_dataset:
                if f not in ['Type','Chromosome']:
                    if f in features_values:
                        if len(features_values[f]) != 0:
                            chrom_cancer_dataset[f].append(stat.median(features_values[f]))
                        else:
                            chrom_cancer_dataset[f].append('')
                    else:
                        chrom_cancer_dataset[f].append('')


    return pd.DataFrame.from_dict(chrom_cancer_dataset)

def add_TSG_OG(Dataset):
    Charm = pd.read_csv("davoli et al (CHARM).csv")
    CHARMscore = {}
    for rowNum, row in Charm.iterrows():
        CHARMscore[row['Chromosome']] = (row['Density_TSG_in_the_arm'], row['Density_OG_in_the_arm'])

    OG_den = []
    TSG_den = []
    for rowNum, row in Dataset.iterrows():
        arm = row['Chromosome']
        TSG_den.append(CHARMscore[arm][0])
        OG_den.append(CHARMscore[arm][1])

    Dataset['TSG_density'] = TSG_den
    Dataset['OG_density'] = OG_den
    return Dataset

def add_labels(Dataset):

    chrom_calls = pd.read_csv("cancer_aneuploidy.csv")

    cancer_calls = {}
    for rowNum, row in chrom_calls.iterrows():
        row = row.to_dict()
        cancer_name = row['Type'].replace('_', ' ').lower()
        if cancer_name not in cancer_calls:
            cancer_calls[cancer_name] = {}

        for item in row:
            if item != 'Type':
                cancer_calls[cancer_name][item] = row[item]
    calls = []

    for rowNum, row in Dataset.iterrows():
        cancer = row['Type'].lower()
        arm = str(int(row['Chromosome']))
        calls.append(cancer_calls[cancer][arm])

    Dataset['Label'] = calls

    return Dataset

def feature_rename(Dataset):
    Dataset = Dataset.rename(columns={'expression': 'Normal tissue expression', 'num interactors': 'Num PPIs', 'num elevated interactors':'Num elevated PPIs', 'num specific interactions':'Num tissue-specific PPIs',
    'diff net med':'Differential PPIs','median tipa pathways': 'Process activity',
    'egene':'Gene eQTL','paralogs ratio highest identity':'Paralogs compensation','fetus development':'Fetus development',
    'childhood development':'Childhood development','young development':'Young development','development cv':'Development variations',
    '(neutarl) TCGA expression':'TCGA expression','(neutral) median crispr score':'CRISPR essentiality score','(neutral) median rnai score':'RNAI essentiality score',
    'tissue specific genes density':'Tissue-specific gene density','essential genes density':'Essential gene density (CRISPR)',
    'OG_density':'Oncogene density','TSG_density':'Tumor suppressor gene density', 'Essentiality density (Nichols)': 'Essential gene density (Nichols)','(neutral) CCL expression':'CCL expression'})
    return Dataset

def construct_Dataset_main():
    print("---------- Constructing cancer type and chromosome arm dataset ----------")

    chromosome_dataset = pd.read_csv('processed_datasets/7_median_top_10_percent.csv')
    gene_dataset = pd.read_csv('processed_datasets/6_Genes z score.csv')
    disease_tissue_association = pd.read_csv("../../Datasets/cancer tissue association.csv")


    features_types = get_features_types(chromosome_dataset)

    chrom_cancer_dataset = assign_fieldnames(features_types)


    chromosome_cancer_dataset = calculation_per_Chrom(chromosome_dataset, disease_tissue_association, chrom_cancer_dataset)  # construct dataset as pairs of cancer type and chromosome arm
    chromosome_cancer_dataset = chromosome_cancer_dataset.drop(columns=[''])
    print('---------- Add TSG and OG densities to dataset ----------')
    chromosome_cancer_dataset = add_TSG_OG(chromosome_cancer_dataset)                                     # add TSG and OG densities features



    print('---------- Calculate tissue-specific gene density and add it to dataset ----------')

    chromosome_cancer_dataset = ts_densities_main(gene_dataset,chromosome_cancer_dataset)         # add tissue-specific gene density feature



    print('---------- Calculate essential gene density based on CRISPR and add it to dataset ----------')

    chromosome_cancer_dataset = crispr_densities_main(gene_dataset,chromosome_cancer_dataset)           # add essential gene density feature


    print('---------- Calculate essential gene density based on Nichols and add it to dataset ----------')
    chromosome_cancer_dataset = add_Essential_gene_density_Nichols(chromosome_cancer_dataset)


    print('---------- Normalize expression features ----------')

    chromosome_cancer_dataset = standardization(chromosome_cancer_dataset)
    chromosome_cancer_dataset['paralogs ratio highest identity'] = chromosome_cancer_dataset['paralogs ratio highest identity']*-1
    chromosome_cancer_dataset = add_labels(chromosome_cancer_dataset)
    chromosome_cancer_dataset = feature_rename(chromosome_cancer_dataset)
    chromosome_cancer_dataset.to_csv('cancer_chrom_dataset.csv',index = False)





def whole_Chromosome_data_construction():
    fix_gene_zscore_file()
    calculate_median_top(pd.read_csv('processed_datasets/6_Genes z score.csv'))
    construct_Dataset_main()

