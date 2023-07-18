import statistics

from calculate_tissue_specific_densities import *
from calculate_essential_densities import *
from essential_gene_density_Nichols import *
from standardize_expression_features import *


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
    fieldnames = {'Arm':[],'Type':[]}

    for f in features_types:
        fieldnames[f] = []
       
    return fieldnames

def get_features_values(row,g_tissue,e_tissue,c_type):
    features_values = {}
    for item in row:
        if item not in ['chromosome arm']:
            item2 = item.split('_')
            
            # Getting current tissue and curr feature
            curr_tissue = item2[0].strip().lower().replace(' ','').replace('-','')
            curr_feature = ''
            for i in range(1,len(item2)):
                curr_feature += item2[i]
                curr_feature += ' ' 
            curr_feature = curr_feature[:-1]
            
            
            if 'tcga' in curr_feature.lower():
                if c_type in curr_tissue:
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

def calculation_per_ChromArm(z_score_dataset, disease_tissue_association, arm_cancer_dataset):

    for rowNum,row in z_score_dataset.iterrows():

        row = row.to_dict()
        if row['chromosome arm'] != '21p':

            for rowNum2, row2 in disease_tissue_association.iterrows():

                arm_cancer_dataset['Arm'].append(row['chromosome arm'])
                arm_cancer_dataset['Type'].append(row2['Type'])

                g_tissue = row2['Tissue - GTEx'].replace(" ","").lower()

                if str(row2['Tissue - DepMap']) != 'nan':
                    e_tissue = row2['Tissue - DepMap'].replace(" ","").lower()
                else:
                    e_tissue = 'nan'

                features_values = get_features_values(row,g_tissue,e_tissue,row2['Type'].lower())



                for f in arm_cancer_dataset:
                    if f not in ['Type','Arm']:
                        if f in features_values:
                            if len(features_values[f]) != 0:
                                arm_cancer_dataset[f].append(statistics.median(features_values[f]))
                            else:
                                arm_cancer_dataset[f].append('')
                        else:
                            arm_cancer_dataset[f].append('')


    return pd.DataFrame.from_dict(arm_cancer_dataset)

def add_labels(Dataset):

    arm_calls = pd.read_csv("../article_results/Figure 1/cancer_arm_aneuploidy.csv")
    cancer_calls = {}
    for rowNum, row in arm_calls.iterrows():
        row = row.to_dict()
        cancer_name = row['Type'].replace('_',' ')
        if cancer_name not in cancer_calls:
            cancer_calls[cancer_name] = {}
        
        for item in row:
            if item!='Type':
                cancer_calls[cancer_name][item]=row[item]
    calls = []

    for rowNum, row in Dataset.iterrows():
        cancer = row['Type']
        arm = row['Arm']
        calls.append(cancer_calls[cancer][arm])
    
    Dataset['Label'] = calls
    
    return Dataset

def add_TSG_OG(Dataset):

    Charm = pd.read_csv("../Datasets/davoli et al (CHARM).csv")
    arm_CHARMscore = {}
    for rowNum, row in Charm.iterrows():
        arm_CHARMscore[row['Arm']] = (row['Density_TSG_in_the_arm'],row['Density_OG_in_the_arm'])

    OG_den = []
    TSG_den = []
    for rowNum, row in Dataset.iterrows():
        arm = row['Arm']
        TSG_den.append(arm_CHARMscore[arm][0])
        OG_den.append(arm_CHARMscore[arm][1])
        
  
    Dataset['TSG_density'] = TSG_den
    Dataset['OG_density'] = OG_den
    return Dataset

def feature_rename(Dataset):
    Dataset = Dataset.rename(columns={'expression': 'Normal tissue expression', 'num interactors': 'Num PPIs', 'num elevated interactors':'Num elevated PPIs', 'num specific interactions':'Num tissue-specific PPIs',
    'diff net med':'Differential PPIs','median tipa pathways': 'Process activity',
    'egene':'Gene eQTL','paralogs ratio highest identity':'Paralogs compensation','fetus development':'Fetus development',
    'childhood development':'Childhood development','young development':'Young development','development cv':'Development variations',
    '(neutarl) tcga expression':'TCGA expression','(neutral) median crispr score':'CRISPR essentiality score','(neutral) median rnai score':'RNAI essentiality score',
    'tissue specific genes density':'Tissue-specific gene density','essential genes density':'Essential gene density (CRISPR)',
    'OG_density':'Oncogene density','TSG_density':'Tumor suppressor gene density', 'Essentiality density (Nichols)': 'Essential gene density (Nichols)'})
    return Dataset

def construct_Dataset_main(Dataset_arm, Dataset_genes_zscore, Dataset_genes):
    print("---------- Constructing cancer type and chromosome arm dataset ----------")


    disease_tissue_association = pd.read_csv("../Datasets/cancer tissue association.csv")
    features_types = get_features_types(Dataset_arm)
    arm_cancer_dataset = assign_fieldnames(features_types)


    arm_cancer_dataset = calculation_per_ChromArm(Dataset_arm, disease_tissue_association, arm_cancer_dataset)  # construct dataset as pairs of cancer type and chromosome arm
    print('---------- Add TSG and OG densities to dataset ----------')
    arm_cancer_dataset = add_TSG_OG(arm_cancer_dataset)                                     # add TSG and OG densities features
    print('---------- Calculate tissue-specific gene density and add it to dataset ----------')
    arm_cancer_dataset = ts_densities_main(Dataset_genes_zscore,arm_cancer_dataset)         # add tissue-specific gene density feature
    print('---------- Calculate essential gene density based on CRISPR and add it to dataset ----------')
    arm_cancer_dataset = crispr_densities_main(Dataset_genes, arm_cancer_dataset)           # add essential gene density feature
    print('---------- Calculate essential gene density based on Nichols and add it to dataset ----------')
    arm_cancer_dataset = add_Essential_gene_density_Nichols(arm_cancer_dataset)
    print('---------- Normalize expression features ----------')
    arm_cancer_dataset = standardization(arm_cancer_dataset)
    arm_cancer_dataset['paralogs ratio highest identity'] = arm_cancer_dataset['paralogs ratio highest identity']*-1
    arm_cancer_dataset = add_labels(arm_cancer_dataset)
    arm_cancer_dataset = feature_rename(arm_cancer_dataset)
    arm_cancer_dataset.to_csv('Processed datasets/8_Arm_Cancer_Dataset.csv', index = False)





