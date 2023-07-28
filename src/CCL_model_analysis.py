import pandas as pd
import statistics
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from sklearn.metrics import precision_recall_curve
from sklearn.metrics import auc
from sklearn.metrics import roc_curve

import xgboost as xgb
from sklearn.ensemble import GradientBoostingClassifier

import shap



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

def null_to_median(Dataset):
    feature_median = {}

    for column in Dataset.columns.tolist():
        values_median = []
        values = Dataset[column].values.tolist()
        for val in values:
            if val == val:
                values_median.append(val)
        feature_median[column] = statistics.median(values_median)
    for rowNum, row in Dataset.iterrows():
        row = row.to_dict()
        for item in row:
            if row[item] != row[item]:
                Dataset.at[rowNum, item] = feature_median[item]


def get_merged_dataframes(train, test):
    merged_dataframe = [train, test]
    merged_dataframe = pd.concat(merged_dataframe)

    indexs = {}
    i = 0
    for ind in list(merged_dataframe.index):
        indexs[ind] = i
        i += 1
    merged_dataframe = merged_dataframe.rename(index=indexs)
    return merged_dataframe

def SHAP_summary_plot(shap_values, X_test, class_model):
    color_map = ["seagreen", "gold", "darkorange"]
    newcmp = LinearSegmentedColormap.from_list('testCmap', colors=color_map, N=256)

    fig = plt.figure(figsize=(10, 8))
    shap.summary_plot(shap_values, X_test, cmap=newcmp, max_display=10, axis_color='black')
    fig.savefig('../article_results/Figure 3/'+class_model+' vs Neutral__ SHAP direction.png', bbox_inches='tight')

def SHAP_bar_plot(shap_values, X_test, class_model):

    fig = plt.figure(figsize=(10, 8))
    shap.summary_plot(shap_values, X_test, max_display=10, plot_type='bar', color='black', axis_color='black')
    fig.savefig('../article_results/Figure 3/'+class_model+' vs Neutral__ SHAP_top_10_features_contribution.png', bbox_inches='tight')

    dict_to_write = {'Feature': [], 'Importance value': []}
    for i in range(0, 16):
        dict_to_write['Feature'].append(np.abs(shap_values).mean(0)[i])
        dict_to_write['Importance value'].append(X_test.columns.tolist()[i])
    pd.DataFrame.from_dict(dict_to_write).to_csv('../article_results/Figure 3/'+class_model+' vs Neutral__ feature importance.csv', index=False)


def SHAP_PLOTS(X, X_test, y,  model,class_model):

    model.fit(X, y)



    # Create object that can calculate shap values
    explainer = shap.TreeExplainer(model)

    # calculate shap values.
    shap_values = explainer.shap_values(X_test)

    SHAP_summary_plot(shap_values, X_test, class_model)
    SHAP_bar_plot(shap_values, X_test, class_model)

def compute_plot_vectors(tests, predictions):
    fpr, tpr, _ = roc_curve(tests, predictions)
    roc_auc = auc(fpr, tpr)

    prec, recall, _ = precision_recall_curve(tests, predictions)
    PR_auc = auc(recall, prec)

    return fpr, tpr, roc_auc, prec, recall, PR_auc

def model_multi_class_tenfold(X_train, y_train, X_test, y_test, model, class_color, threshold, class_model):

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 8))
    lw = 2
    ax1.set_xlim([0.0, 1.0])
    ax1.set_ylim([0.0, 1.05])

    ax2.set_xlim([0.0, 1.0])
    ax2.set_ylim([0.0, 1.05])

    # We will do leave one out for the test and train the model on the rest.


    predictions = model.fit(X_train, y_train).predict_proba(X_test)[:, 1]
    fpr, tpr, roc_auc, prec, recall, PR_auc = compute_plot_vectors(np.asarray(y_test.tolist()), predictions)
    ax1.plot(fpr, tpr, lw=lw, color=class_color)
    ax2.plot(recall, prec, lw=lw, color=class_color)

    ax1.plot([0, 1], [0, 1], color='black', lw=lw, linestyle='--', label = 'auc = %.2f '%roc_auc)
    ax1.set_xlabel('False Positive Rate')
    ax1.set_ylabel('True Positive Rate')
    ax1.set_title('ROC curve')
    ax1.legend(loc="lower right")
    ax2.plot([0, 1], [threshold, threshold], color='black', lw=lw, linestyle='--', label = 'prc = %.2f '%PR_auc)
    ax2.set_xlabel('Recall')
    ax2.set_ylabel('Percision')
    ax2.set_title('Percision-Recall curve')
    ax2.legend(loc="lower right")

    plt.savefig('../article_results/Figure 3/Panel A - model performance_'+class_model+'.png', bbox_inches='tight')

def drop_row_based_on_class(dataset,class_to_drop):
    drop_rows = []
    for rowNum, row in dataset.iterrows():
        if row['Label'] == class_to_drop:
            drop_rows.append(rowNum)
    return drop_rows

def calculate_CCL_models():


    ### Importing and preprocessing TCGA dataset
    TCGA_dataset = pd.read_csv("Processed datasets/8_Arm_Cancer_Dataset.csv")
    columns_to_drop = ['Development variations', 'Childhood development', 'Fetus development', 'Young development']
    TCGA_dataset = TCGA_dataset.drop(columns=columns_to_drop)
    TCGA_dataset = TCGA_dataset.rename(columns={'TCGA expression':'CCL expression'})

    ### Importing and preprocessing CCLE dataset
    CCL_dataset = pd.read_csv("Processed datasets/9_CCL_arm_aneuploidy.csv")

    ### Picking to Classes to compare
    first_class = ['Gain','Loss']
    second_class = 'Neutral'


    for c1 in first_class:
        label_list = list(TCGA_dataset['Label']) + list(CCL_dataset['Label'])
        counter = {x: label_list.count(x) for x in label_list}
        counter['total'] = len(label_list)

        labels, class_to_drop, threshold, total_count = get_data_info(c1, second_class, counter)

        ### Preprocessing the required dataset regarding to the classes that were picked
        TCGA_rows_to_drop = drop_row_based_on_class(TCGA_dataset,class_to_drop)
        TCGA_class_dataset = TCGA_dataset.drop(TCGA_rows_to_drop)
        CCL_rows_to_drop = drop_row_based_on_class(CCL_dataset,class_to_drop)
        CCL_class_dataset = CCL_dataset.drop(CCL_rows_to_drop)


        columns_to_drop = ['Type', 'Arm']

        TCGA_class_dataset = TCGA_class_dataset.drop(columns=columns_to_drop)
        TCGA_class_dataset= TCGA_class_dataset.replace({'Label': labels})


        CCL_class_dataset = CCL_class_dataset.drop(columns=columns_to_drop)
        CCL_class_dataset = CCL_class_dataset.replace({'Label': labels})

        total_count = len(CCL_class_dataset['Label'].tolist())
        CCL_class_dataset = CCL_class_dataset.set_axis(range(0, total_count))

        ###Imputation
        null_to_median(TCGA_class_dataset)
        null_to_median(CCL_class_dataset)

        ### train and test datasets
        X_train = TCGA_class_dataset.loc[:, TCGA_class_dataset.columns != 'Label']
        y_train = TCGA_class_dataset['Label']

        X_test = CCL_class_dataset.loc[:, CCL_class_dataset.columns != 'Label']
        y_test = CCL_class_dataset['Label']
        X_test.insert(8, 'CCL expression', X_test.pop('CCL expression'))
        X = get_merged_dataframes(X_train, X_test)
        y = get_merged_dataframes(y_train, y_test)

        ### Pick model and Run SHAP analysis
        if c1 == 'Loss':
            model = xgb.XGBClassifier(eval_metric='mlogloss')
            class_color = '#4472C4'
            threshold = 0.11

        if c1 == 'Gain':
            model = GradientBoostingClassifier()
            class_color = '#BF0000'
            threshold = 0.15
        SHAP_PLOTS(X, X_test, y,model,c1)
        model_multi_class_tenfold(X_train, y_train, X_test, y_test, model, class_color, threshold, c1)
