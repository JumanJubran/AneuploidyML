import pandas as pd

def find_range_val_per_cancer(feature,final_dataset):
    max_per_cancer = {}
    min_per_cancer = {}
    for rowNum, row in final_dataset.iterrows():
        if row['Type'] not in max_per_cancer:
            max_per_cancer[row['Type']] = 0
            min_per_cancer[row['Type']] = 100
        
        if row[feature] > max_per_cancer[row['Type']]:
            max_per_cancer[row['Type']] = row[feature]
        
        if row[feature] < min_per_cancer[row['Type']]:
            min_per_cancer[row['Type']] = row[feature]
    return max_per_cancer, min_per_cancer


def standardization(Dataset):

    features_to_standardize = ['expression','(neutarl) TCGA expression']

    for feature in features_to_standardize:
        stand_exp = []
        max_per_cancer, min_per_cancer = find_range_val_per_cancer(feature,Dataset)
        for rowNum, row in Dataset.iterrows():
            stand_exp.append((row[feature]-min_per_cancer[row['Type']])/(max_per_cancer[row['Type']]-min_per_cancer[row['Type']]))
        #f_name = feature+'_standardized'
        Dataset[feature] = stand_exp
    return Dataset
