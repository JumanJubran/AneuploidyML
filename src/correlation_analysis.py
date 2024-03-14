import statistics
from get_info import *

import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import pandas as pd
import numpy as np



def get_arm_frequencies(data_type):
    cancer_types = get_cancer_types_in_model(data_type)
    freq_per_arm = {}

    for t in cancer_types:
        if data_type == 'TCGA':
            freq_tables = pd.read_csv("../Datasets/GISTIC_aneuploidy_stats_per_cancer/"+t+".csv")
            del_freq = 'Del Frequency'
            amp_freq = 'Amp Frequency'
        else:
            freq_tables = pd.read_csv("../Datasets/CCL data/statistics per cancer/"+t+".csv")
            del_freq = 'Del frequency - WGD'
            amp_freq = 'Amp frequency - WGD'

        for rowNum, row in freq_tables.iterrows():

            arm = row['Arm']
            if arm == arm:
                if arm not in freq_per_arm:
                    freq_per_arm[arm] = [[], []]

                freq_per_arm[arm][0].append(row[del_freq])
                freq_per_arm[arm][1].append(row[amp_freq])
    return freq_per_arm

def get_arm_features_values(dataset):
    feature_per_arm = {}
    for rowNum, row in dataset.iterrows():
        row = row.to_dict()

        if row['Arm'] not in feature_per_arm:
            feature_per_arm[row['Arm']] = {}
            flag = 0
        else:
            flag = 1
        for feature in row:
            if feature not in ['Arm', 'Type', 'Label']:
                if flag == 0:
                    feature_per_arm[row['Arm']][feature] = []
                if row[feature] == row[feature]:
                    feature_per_arm[row['Arm']][feature].append(row[feature])
    return feature_per_arm

def calculate_median_per_arm(dataset, data_type):

    freq_per_arm = get_arm_frequencies(data_type)
    feature_per_arm = get_arm_features_values(dataset)

    dict_to_write = {'Arm': [], 'Del frequency': [], 'Amp frequency': []}

    for arm in freq_per_arm:
        dict_to_write['Arm'].append(arm)
        dict_to_write['Del frequency'].append(statistics.median(freq_per_arm[arm][0]))
        dict_to_write['Amp frequency'].append(statistics.median(freq_per_arm[arm][1]))
        features = feature_per_arm[arm]
        for feature in features:
            if feature not in dict_to_write:
                dict_to_write[feature] = []
            dict_to_write[feature].append(statistics.median(features[feature]))

    data_to_write = pd.DataFrame.from_dict(dict_to_write)
    data_to_write.to_csv( "Processed datasets/"+data_type+"_median_per_arm.csv",index=False)
    return data_to_write


def add_row_to_file(dict_to_write, freq_n, feature_n, m, b, r, p):
    dict_to_write['Frequency type'].append(freq_n)
    dict_to_write['Feature'].append(feature_n)
    dict_to_write['m'].append(m)
    dict_to_write['b'].append(b)
    dict_to_write['r'].append(r)
    dict_to_write['p-value'].append(p)
    return dict_to_write

def scatterplot(dataset, feature_n, freq_n, x, y,color,data_type):
    sns.set(rc={'figure.figsize': (18, 12)})
    sns.set(style='ticks', font_scale=1.2)

    plt.xlim((min(x) - ((max(x) - min(x)) / 100 * 4), max(x) + ((max(x) - min(x)) / 100 * 4)))

    g = sns.regplot(x=x, y=y, data=dataset, scatter_kws={"color": "#67879B"}, line_kws={"color": color})

    plt.xlabel(feature_n, fontsize=28)
    plt.ylabel(freq_n, fontsize=28)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)

    ### Add arm label of each dot
    for line in range(0, dataset.shape[0]):
        g.text(dataset[feature_n][line], dataset[freq_n][line] + 0.0001, dataset['Arm'][line],
               horizontalalignment='left', fontsize=20, color='black')

    ### Save plot
    if feature_n in ['Tumor suppressor gene density','Oncogene density','TCGA expression','CCL expression','Essential gene density (Nichols)']:
        if data_type == 'TCGA':
            plt.savefig("../article_results/Figure 2/Panel E_" + freq_n + "_" + feature_n + "_TCGA_correlation.png")
        else:
            plt.savefig("../article_results/Figure 3/Panel F_" + freq_n + "_" + feature_n + "_CCL_correlation.png")

    elif feature_n == 'Paralogs compensation':
        if freq_n == 'Del frequency':
            plt.savefig("../article_results/Figure 5/Panel A_" + freq_n + "_" + feature_n + "_"+data_type+"_correlation.png")

    else:
        if data_type == 'TCGA':
            plt.savefig( "../article_results/Figure S6/" + freq_n + "_" + feature_n + "_" + data_type + "_correlation.png")
    plt.clf()


def correlation(dataset,data_type):
    low_color = '#2E8B57'
    high_color = '#ff8c00'

    freq_type = ['Del frequency', 'Amp frequency']

    dict_to_write = {'Frequency type': [], 'Feature': [], 'm': [], 'b': [], 'r': [], 'p-value': []}

    for feature_n in dataset.columns.tolist():
        if feature_n not in ['Arm', 'Del frequency', 'Amp frequency']:
            for freq_n in freq_type:

                y = dataset[freq_n].tolist()  ### Frequency's values
                x = dataset[feature_n].tolist()  ### Feature's values

                m, b = np.polyfit(np.asarray(x), np.asarray(y), 1)
                r, p = spearmanr(y, x)

                if abs(r) < 0.1:
                    color = '#808080'
                elif r < 0:
                    color = low_color
                else:
                    color = high_color

                ### Scatterplot
                scatterplot(dataset, feature_n, freq_n, x, y,color,data_type)

                ### Add the result of the specific correlation to the dictionary
                dict_to_write = add_row_to_file(dict_to_write, freq_n, feature_n, m, b, r, p)

    ### Save the correlation results
    if data_type == 'TCGA':
        pd.DataFrame.from_dict(dict_to_write).to_csv("../article_results/Figure 2/correlation_results.csv", index=False)
    else:
        pd.DataFrame.from_dict(dict_to_write).to_csv("../article_results/Figure 3/correlation_results.csv", index=False)


def calculate_correlations():
    median_per_arm = calculate_median_per_arm(pd.read_csv("Processed datasets/8_Arm_Cancer_Dataset.csv"), 'TCGA')
    correlation(median_per_arm,'TCGA')
    median_per_arm = calculate_median_per_arm(pd.read_csv("Processed datasets/9_CCL_arm_aneuploidy.csv"), 'CCL')
    correlation(median_per_arm, 'CCL')
