import statistics as stat

from get_info import *


def z_score_per_gene(features_values,gene):
    gene_z_score = {}
    for feature in features_values:
        

        gene_values = list(features_values[feature].values())

        gene_values = [x for x in gene_values if str(x) != 'nan' and x != '']

        if len(gene_values)>1:
            gene_z_score[feature] = {}

            mean = stat.mean(gene_values)
            stdv = stat.stdev(gene_values)

            for tissue in features_values[feature]:

                if stdv == 0:
                    if features_values[feature][tissue] == features_values[feature][tissue]:

                        gene_z_score[feature][tissue] = 0
                    else:
                        gene_z_score[feature][tissue] = ''
                else:
                    gene_score = features_values[feature][tissue]
                    if gene_score != '':

                        z_score = (gene_score - mean)/stdv
                        gene_z_score[feature][tissue] = z_score
                    else:
                        gene_z_score[feature][tissue] = ''


        else:
            gene_z_score[feature] = {}
            for tissue in features_values[feature]:
                gene_z_score[feature][tissue] = ''
                

    return gene_z_score

def add_zscore_to_dict(z_score_dict,gene_z_score):
    for feature in gene_z_score:
        for tissue in gene_z_score[feature]:
            f_t = tissue+'_'+feature
            z_score_dict[f_t].append(gene_z_score[feature][tissue])
    return(z_score_dict)

def essential_opposite_direction(Dataset):
    '''
    Changing the directionality of crispr score so the high values indicate essentiality of a gene
    :param Dataset:
    :return:
    '''
    for item in Dataset.columns.tolist():
        if 'crispr' in item or 'rnai' in item:
            Dataset[item] = [element * (-1) for element in Dataset[item]]
    return Dataset

def calculate_z_score(Dataset):

    z_score_dict = {}

    for item in Dataset.columns.tolist():
        z_score_dict[item] = []
       
    for rowNum, row in Dataset.iterrows():
        features_values = {}
        row = row.to_dict()


        for item in row:
            if item not in ['Chromosome_arm','Gene_ID','Gene_name']:
                if 'diff' not in item and 'dif' not in item and 'tipa' not in item and 'preferential' not in item:

                    tissue_name, feature_name = split_name(item)

                    if feature_name not in features_values:
                        features_values[feature_name] = {}

                    features_values[feature_name][tissue_name] = row[item]

                else:
                    z_score_dict[item].append(row[item])
            else:
                z_score_dict[item].append(row[item])
     
  
        gene_z_score = z_score_per_gene(features_values, row['Gene_name'])
        
        z_score_dict = add_zscore_to_dict(z_score_dict, gene_z_score)
        

    Dataset = pd.DataFrame.from_dict(z_score_dict)
    Dataset = essential_opposite_direction(Dataset)


    return Dataset
