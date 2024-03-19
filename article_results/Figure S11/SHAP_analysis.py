import pandas as pd
import statistics
import numpy as np

### ML Models and SHAP
import xgboost as xgb
import shap
from sklearn.ensemble import RandomForestClassifier

### matplotlib functions
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

'''
Getting data Info
'''
def get_data_info(Class, second_class, counter):
    total_count = counter['total']
    if Class == 'Loss':
        class_to_drop = 1
        labels = {-1: 1, 0: 0, 1: 0}
        threshold = counter[-1] / counter['total']
        if second_class == 'Neutral':
            total_count = counter[-1] + counter[0]
            threshold = counter[-1] / total_count
    else:
        class_to_drop = -1
        labels = {-1: 0, 0: 0, 1: 1}
        threshold = counter[1] / counter['total']
        if second_class == 'Neutral':
            total_count = counter[1] + counter[0]
            threshold = counter[1] / total_count

    return labels, class_to_drop, threshold, total_count


'''
Imputation: To median
'''
def null_to_median(dataset_per_model):
    feature_median = {}

    for column in dataset_per_model.columns.tolist():
        values_median = []
        values = dataset_per_model[column].values.tolist()
        for val in values:
            if val == val:
                values_median.append(val)
        feature_median[column] = statistics.median(values_median)
    for rowNum, row in dataset_per_model.iterrows():
        row = row.to_dict()
        for item in row:
            if row[item] != row[item]:
                dataset_per_model.at[rowNum, item] = feature_median[item]

    return dataset_per_model


def SHAP_summary_plot(shap_values, X_test, file_name):
    color_map = ["seagreen", "gold", "darkorange"]
    newcmp = LinearSegmentedColormap.from_list('testCmap', colors=color_map, N=256)

    fig = plt.figure(figsize=(10, 8))
    shap.summary_plot(shap_values, X_test, cmap=newcmp, max_display=10, axis_color='black')


    fig.savefig('results/'+file_name + '__SHAP_direction.png', bbox_inches='tight')



def SHAP_bar_plot(shap_values, X_test, file_name):
    fig = plt.figure(figsize=(10, 8))
    shap.summary_plot(shap_values, X_test,plot_type='bar', color='black', axis_color='black')


    fig.savefig('results/'+file_name + '__SHAP_all_features_contribution.png', bbox_inches='tight')

    fig = plt.figure(figsize=(10, 8))
    shap.summary_plot(shap_values, X_test, max_display=10, plot_type='bar', color='black', axis_color='black')

    fig.savefig('results/'+file_name + '__SHAP_top10_features_contribution.png',bbox_inches='tight')



def save_features_importance_as_output_file(X, shap_values, file_name,instance_labels):
    dict_to_write = {'Feature': [], 'Importance value': []}
    shap_values_dict = {'Instance': instance_labels}
    for i in range(0, 20):

        dict_to_write['Importance value'].append(np.abs(shap_values).mean(0)[i])
        dict_to_write['Feature'].append(X.columns.tolist()[i])


        shap_values_dict[X.columns.tolist()[i]] = shap_values[:, i]

        pd.DataFrame.from_dict(shap_values_dict).to_csv('results/'+file_name + '__feature_shap_values.csv', index=False)
        pd.DataFrame.from_dict(dict_to_write).to_csv('results/'+file_name + '__feature_importance_valuess.csv',index=False)

'''
SHAP analysis of top contributing features
'''
def SHAP_PLOTS(X, y, model, file_name,instance_labels):
    model.fit(X, y)

    # Create object that can calculate shap values
    explainer = shap.TreeExplainer(model)

    # calculate shap values. Which showes the direction and contribution of each feature

    shap_values = explainer.shap_values(X)

    if 'Gain' in file_name:
        shap_values = shap_values[1]

    save_features_importance_as_output_file(X, shap_values, file_name,instance_labels)
    SHAP_summary_plot(shap_values, X, file_name)
    SHAP_bar_plot(shap_values, X, file_name)


'''
process the data based on the model. There are 2 models, including: 
Gain vs. Neutral, 
Loss vs. Neutral, 
'''
def process_dataset_per_model(data, c1, c2):
    num_observation_per_class = {x: list(data['Label']).count(x) for x in list(data['Label'])}
    num_observation_per_class['total'] = len(list(data['Label']))

    labels, class_to_drop, threshold, total_count = get_data_info(c1, c2, num_observation_per_class)

    if c2 == 'Neutral':
        for rowNum, row in data.iterrows():
            if row['Label'] == class_to_drop:
                data = data.drop(rowNum)

    data = data.replace({'Label': labels})
    data = data.set_axis(range(0, total_count))

    return data, threshold


'''
Main Function
'''
def SHAP_analysis_of_features_contributing():
    Dataset = pd.read_csv("cancer_chrom_dataset.csv")

    Class = ['Gain', 'Loss']

    best_performing_model_per_model = {('Gain', 'Neutral'): RandomForestClassifier(),
                                       ('Loss', 'Neutral'): xgb.XGBClassifier(use_label_encoder=False, eval_metric='logloss')}

    ### Choose classes to compare

    for c1 in Class:
        c2 = 'Neutral'
        ### Directory of the output
        file_name = c1 + ' vs ' + c2



        dataset_per_model, threshold = process_dataset_per_model(Dataset, c1, c2)

        instance_labels = []
        for rowNum, row in dataset_per_model.iterrows():
            instance_labels.append(row['Type'] + '-' + (str(row['Chromosome'])))

        dataset_per_model = dataset_per_model.drop(columns=['Type', 'Chromosome'])

        dataset_per_model = null_to_median(dataset_per_model)

        X = dataset_per_model.loc[:, dataset_per_model.columns != 'Label']
        y = dataset_per_model['Label']

        model = best_performing_model_per_model[(c1, c2)]

        SHAP_PLOTS(X, y, model, file_name,instance_labels)

