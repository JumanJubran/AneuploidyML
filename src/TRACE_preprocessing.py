import pandas as pd
import statistics as stat

def remove_embedding_features_from_TRACE(Trace_Dataset):
    '''
    Description: filtering embedding features from TRACE.
    param: Trace_Dataset:
    return: Trace Dataset without embeddings
    '''


    wanted_features = []
    unwanted_features = []
    for i in Trace_Dataset.columns.tolist():
        if 'embedding' not in i and 'HP' not in i:
            wanted_features.append(i)
        else:
            unwanted_features.append(i)

    Trace_Dataset_without_embeddings = Trace_Dataset.drop(columns=unwanted_features)
    return Trace_Dataset_without_embeddings

def add_chromosomes_arm(Trace_Dataset):
    '''
    Description: adding chromosome arm regarding genes tp TRACE datset
    param: Trace_Dataset
    return: Trace Dataset with chromosome arm
    '''

    protein_per_chromosome = pd.read_csv("../Datasets/protein coding genes - chromosomal location(bioMart).csv")

    proteins_chromosome = {}

    for rowNum,row in protein_per_chromosome.iterrows():
        chromosome = str(row['Chromosome/scaffold name'])+str(row['Karyotype band'])[0]
        proteins_chromosome[row['Gene stable ID']] = (row['Gene name'],chromosome)



    gene_name = []
    chromosome_arm = []
    drop_rows = []
    for rowNum, row in Trace_Dataset.iterrows():

        if row['Unnamed: 0'] in proteins_chromosome:
            chromosome_arm.append(proteins_chromosome[row['Unnamed: 0']][1])
            gene_name.append(proteins_chromosome[row['Unnamed: 0']][0])
        else:

            chromosome_arm.append('')
            gene_name.append('')
            drop_rows.append(rowNum)

    Trace_Dataset.insert(1, 'Gene_name', gene_name)
    Trace_Dataset.insert(2, 'Chromosome_arm', chromosome_arm)
    Trace_Dataset = Trace_Dataset.drop(drop_rows)
    Trace_Dataset = Trace_Dataset.set_axis(range(0, len(list(Trace_Dataset['Gene_name']))))

    return Trace_Dataset

def fix_features_names(feature_name):
    '''
    Description: fixing the name of the features so it will be easier to extract the tissue and the feature
    param: feature_name
    return: fixed feature name
    '''
    feature_name = feature_name.lower()
    feature_name = feature_name.replace('_-_', '_') \
        .replace('.', '_') \
        .replace(')', '') \
        .replace('(', '') \
        .replace('whole_brain', 'whole-brain') \
        .replace('whole_blood', 'whole-blood') \
        .replace('skin_sun_exposed_lower_leg', 'skin-sun-exposed-lower-leg') \
        .replace('skin_not_sun_exposed_suprapubic', 'skin-not-sun-exposed-suprapubic') \
        .replace('minor_salivary_gland', 'minor-salivary-gland') \
        .replace('heart_left_ventricle', 'heart-left-ventricle') \
        .replace('heart_atrial_appendage', 'heart-atrial-appendage') \
        .replace('esophagus_gastroesophageal_junction', 'esophagus-gastroesophageal-junction') \
        .replace('esophagus_muscularis', 'esophagus-muscularis') \
        .replace('esophagus_mucosa', 'esophagus-mucosa') \
        .replace('cells_transformed_fibroblasts', 'cells-transformed-fibroblasts') \
        .replace('cells_ebv-transformed_lymphocytes', 'cells-ebv-transformed-lymphocytes') \
        .replace('breast_mammary_tissue', 'breast-mammary-tissue') \
        .replace('brain_substantia_nigra', 'brain-substantia-nigra') \
        .replace('brain_spinal_cord_cervical_c-1', 'brain-spinal-cord-cervical-c-1') \
        .replace('brain_putamen_basal_ganglia', 'brain-putamen-basal-ganglia') \
        .replace('brain_nucleus_accumbens_basal_ganglia', 'brain-nucleus-accumbens-basal-ganglia') \
        .replace('brain_frontal_cortex_ba9', 'brain-frontal-cortex-ba9') \
        .replace('brain_cerebellar_hemisphere', 'brain-cerebellar-hemisphere') \
        .replace('brain_caudate_basal_ganglia', 'brain-caudate-basal-ganglia') \
        .replace('brain_anterior_cingulate_cortex_ba24', 'brain-anterior-cingulate-cortex-ba24') \
        .replace('adipose_visceral_omentum', 'adipose-visceral-omentum') \
        .replace('adipose_subcutaneous', 'adipose-subcutaneous') \
        .replace('adipose_subcutaneous', 'adipose-subcutaneous') \
        .replace('muscle_skeletal', 'muscle-skeletal') \
        .replace('nerve_tibial', 'nerve-tibial') \
        .replace('colon_sigmoid', 'colon-sigmoid') \
        .replace('colon_transverse', 'colon-transverse') \
        .replace('brain_hypothalamus', 'brain-hypothalamus') \
        .replace('brain_hippocampus', 'brain-hippocampus') \
        .replace('brain_cortex', 'brain-cortex') \
        .replace('brain_cerebellum', 'brain-cerebellum') \
        .replace('brain_amygdala', 'brain-amygdala') \
        .replace('artery_coronary', 'artery-coronary') \
        .replace('artery_aorta', 'artery-aorta') \
        .replace('artery_tibial', 'artery-tibial') \
        .replace('adrenal_gland', 'adrenal-gland') \
        .replace('small_intestine_terminal_ileum', 'small-intestine-terminal-ileum') \
        .replace('brain-not_specific', 'brain-not-specific') \
        .replace('basal_ganglia', 'basal-ganglia') \
        .replace('spinal_cord', 'spinal-cord') \
        .replace('brain_other', 'brain-other') \
        .replace('whole-brain_lower', 'whole-brain-lower')
    return feature_name

def drop_unwanted_features(Trace_Dataset):
    '''
    Description: droping unwanted features from Trace dataset
    param: Trace_Dataset
    return: Trace dataset with specific features
    '''
    features_to_drop = ['preferential_expression',
                       'num_interactors_dif_med', 'num_interactors_dif_mean',
                       'num_elevated_interactors_dif_median', 'num_elevated_interactors_dif_mean',
                       'num_specific_interactions_dif_median', 'num_specific_interactions_dif_mean',
                       'diff_net_max', 'diff_net_min', 'diff_net_mean',
                       'max_tipa_pathways', 'min_tipa_pathways', 'mean_tipa_pathways',
                       'lcv', 'paralogs_ratio_all', 'causal',
                       'youngteenager_development',
                       'oldteenager_development', 'paths_num_pathways','senior_development','youngmidage_development','oldermidage_development']

    drop_features = []
    fixed_feature_names = {'Unnamed: 0':'Gene_ID','Gene_name':'Gene_name','Chromosome_arm':'Chromosome_arm'}
    for feature_name in Trace_Dataset.columns.tolist():
        if feature_name not in fixed_feature_names:
            fixed_feature_name = fix_features_names(feature_name)

            feature_name_s = fixed_feature_name.split('_')
            feature = ''
            for i in range(1, len(feature_name_s)):
                feature += feature_name_s[i]
                feature += '_'
            feature = feature[:-1]
            if feature in features_to_drop:
                drop_features.append(feature_name)
            else:
                fixed_feature_names[feature_name] = fixed_feature_name

    Trace_Dataset = Trace_Dataset.drop(columns=drop_features)
    Trace_Dataset = Trace_Dataset.rename(columns=fixed_feature_names)
    return Trace_Dataset

def features_development_group(Trace_Dataset):
    '''
    Description: Dividing development features to 3 groups including fetus, childhood, and young.
    param: Trace_Dataset
    return: features of each group
    '''
    fetus_features = []
    childhood_features = ['newborn_development', 'infant_development', 'toddler_development']
    young_features = ['school_development', 'teenager_development','youngteenager_development','oldteenager_development',  'youngadult_development']

    for item in Trace_Dataset.columns.tolist():
        if 'wpc' in item.lower():
            item = item.lower().split('_')
            feature = ''
            for i in range(1, len(item)):
                feature += item[i]
                feature += '_'
            if feature[:-1] not in fetus_features:
                fetus_features.append(feature[:-1])

    return fetus_features, childhood_features, young_features

def get_feature_tissue_pair(column_name):
    '''
    Description: split the column name of the dataset to tissue name and feature type
    param: item
    return: tissue and feature of the related column
    '''
    column_name_s = column_name.split('_')
    feature = ''
    for i in range(1,len(column_name_s)):
        feature += column_name_s[i]
        feature += '_'
    feature = feature[:-1]
    tissue = column_name_s[0]
    return feature, tissue

def add_values(feature, tissue, features_list, development_dict, value):
    '''
    Description: adding the value of the tissue expression to the specific development group
    param: feature, tissue, features_list, development_dict, value:
    return:
    '''
    if feature in features_list:
        if tissue not in development_dict:
            development_dict[tissue] = []
        if value == value:
            development_dict[tissue].append(value)
    return development_dict

def calculate_tissue_median(development_dict, median_dict):
    '''
    Description: calculating the tissue median expression of the development group
    param development_dict, median_dict
    :return: median expression of development groups per tissue
    '''
    for tissue in development_dict:
        if tissue not in median_dict:
            median_dict[tissue] = []
        if len(development_dict[tissue]) != 0:
            median_dict[tissue].append(stat.median(development_dict[tissue]))
        else:
            median_dict[tissue].append('')
    return median_dict

def remove_development_features(Trace_dataset):
    '''
    escription: removing unwanted development features of TRACE
    param: Trace_dataset
    return: Trace dataset without development sub-group features
    '''
    cols_to_drop = []
    for item in Trace_dataset.columns.tolist():
        item_l = item.lower()
        if 'development' in item_l:
            if 'cv' not in item_l:
                cols_to_drop.append(item)
    Trace_dataset = Trace_dataset.drop(columns = cols_to_drop)
    return Trace_dataset

def add_new_development_features(Trace_dataset, median_dict, feature_name):
    '''
    Description: adding the new integrated development features
    param Trace_dataset, median_dict, feature_name
    return: Trace dataset with the intregrated development feature
    '''
    for tissue in median_dict:
        Trace_dataset[tissue+feature_name] = median_dict[tissue]
    return Trace_dataset

def unit_development_features(Trace_Dataset):
    '''
    Description: uniting development features into three categories. Each gene will get the median value the sub-categories expression
    param: Trace_Dataset
    return: Trace_Datsaet with integrated 3 development expression categories
    '''

    # These dicts will hold tissue median expression of the development group (dict[tissue]=[list of median expression of all genes])
    fetus_median_vals = {}
    childhood_median_vals = {}
    young_median_vals = {}


    fetus_features, childhood_features, young_features = features_development_group(Trace_Dataset)

    for rowNum, row in Trace_Dataset.iterrows():
        # These dicts will hold the tissue expressions of the development group of the current gene
        fetus_development = {}
        childhood_development = {}
        young_development = {}

        row = row.to_dict()
        for item in row:
            item_l = item.lower()
            if 'development' in item_l:
                feature, tissue = get_feature_tissue_pair(item_l)

                # Add the expression value to the correct development group
                fetus_development = add_values(feature, tissue, fetus_features, fetus_development, row[item])
                childhood_development = add_values(feature, tissue, childhood_features, childhood_development,row[item])
                young_development = add_values(feature, tissue, young_features, young_development, row[item])


        # Calculate the median
        fetus_median_vals = calculate_tissue_median(fetus_development, fetus_median_vals)
        childhood_median_vals = calculate_tissue_median(childhood_development, childhood_median_vals)
        young_median_vals = calculate_tissue_median(young_development, young_median_vals)


    Trace_Dataset = remove_development_features(Trace_Dataset)

    Trace_Dataset = add_new_development_features(Trace_Dataset, fetus_median_vals, '_fetus_development')
    Trace_Dataset = add_new_development_features(Trace_Dataset, childhood_median_vals, '_childhood_development')
    Trace_Dataset = add_new_development_features(Trace_Dataset, young_median_vals, '_young_development')

    return Trace_Dataset


'''
    ### MAIN FUNCTION! ###
''' 
def TRACE_preprocessing_Main():
    Trace_Dataset = pd.read_csv("../Datasets/TRACE_dataset.csv")

    # Remove embedding features from TRACE dataset
    Trace_Dataset = remove_embedding_features_from_TRACE(Trace_Dataset)

    # Adding chromosome-arm of the genes
    Trace_Dataset = add_chromosomes_arm(Trace_Dataset)

    # Removing unwanted features
    Trace_Dataset = drop_unwanted_features(Trace_Dataset)
    Trace_Dataset.to_csv('Processed datasets/1_TRACE_Dataset_processed_features.csv', index = False)

    # Uniting development features
    Trace_Dataset = unit_development_features(Trace_Dataset)
    Trace_Dataset.to_csv('Processed datasets/2_TRACE_Dataset_united_development_features.csv', index=False)
    return Trace_Dataset