import pandas as pd
from SHAP_analysis import *

def top_driver_pairs(data, label):

    ### 5 topmost features of gain model
    if label == 'Gain':
        top_features = ['Instance', 'Tumor suppressor gene density', 'Oncogene density',
                        'Essential gene density (Nichols)', 'Num PPIs', 'TCGA expression']

    ### 5 topmost features of loss model
    else:
        top_features = ['Instance', 'Tumor suppressor gene density', 'Oncogene density', 'Paralogs compensation',
                        'Num elevated PPIs', 'TCGA expression']


    dict_to_write = {'Pair': [], 'Count (All features)': [], 'Count (Top-five)': []}
    features_per_pair = {}

    for item in list(data['Instance']):
        features_per_pair[item] = []

    ### count how many times the feature's value was among the 10% highest values and 10% lowest values
    for feature in data.columns:
        if feature != 'Instance':
            dict_to_write[feature] = []
            feature_dict = data.set_index('Instance')[feature].to_dict()
            sorted_feature_dict = sorted(feature_dict.items(), key=lambda x: x[1])
            size_of_list = len(sorted_feature_dict)

            ten_per_size = int(size_of_list * 0.1)
            neg_ten_per_size = -1 * ten_per_size
            bottom_ten = sorted_feature_dict[:ten_per_size]

            top_ten = sorted_feature_dict[neg_ten_per_size:]
            extreme_val = bottom_ten + top_ten

            for item in extreme_val:
                features_per_pair[item[0]].append(feature)


    for item in features_per_pair:
        count_all = 0
        count_top_5 = 0
        dict_to_write['Pair'].append(item)
        for feature in dict_to_write:
            if feature not in ['Pair', 'Count (All features)', 'Count (Top-five)']:
                if feature in features_per_pair[item]:
                    dict_to_write[feature].append(1)
                    if feature in top_features:
                        count_top_5 += 1
                    count_all += 1
                else:
                    dict_to_write[feature].append(0)

        dict_to_write['Count (All features)'].append(count_all)
        dict_to_write['Count (Top-five)'].append(count_top_5)

    pd.DataFrame.from_dict(dict_to_write).to_csv(label + " vs. Neutral extreme pairs.csv",index=False)


for label in ['Loss', 'Gain']:
    data = pd.read_csv('../Figure 2/' + label + ' vs Neutral__feature_shap_values.csv')
    top_driver_pairs(data, label)

SHAP_analysis_of_features_contributing()