import csv
import matplotlib.pyplot as plt
import seaborn as sns

from get_info import *



def get_aneuploidy_status(amp_q_val, del_q_val, amp_freq, del_freq, threshold):
    if amp_q_val < threshold and del_q_val < threshold:
        if amp_q_val == del_q_val:
            if amp_freq > del_freq:
                return 1
            else:
                return -1
        elif amp_q_val < del_q_val:
            return 1
        else:
            return -1

    elif amp_q_val < threshold:
        return 1
    elif del_q_val < threshold:
        return -1
    else:
        return 0

def get_qvalue_threshold(data_type):
    if data_type == 'TCGA':
        return 0.05
    elif data_type == 'CCL':
        return 0.15

def results_as_heatmap():
    data = pd.read_csv("../article_results/Figure 1/Panel A _ cancer_arm_aneuploidy.csv")

    data = data.set_index(data['Type'])
    data = data.drop(columns=['Type'])
    fig, ax = plt.subplots(figsize=(12, 8))

    ax.plot()

    g = sns.heatmap(data, cmap='vlag', square=True, linewidth=0.01, linecolor="#222")
    g.set_xticklabels(g.get_xmajorticklabels(), fontsize=14)
    g.set_yticklabels(g.get_ymajorticklabels(), fontsize=16)

    plt.savefig("../article_results/Figure 1/Panel A _ cancer_arm_aneuploidy_heatmap.png",
                bbox_inches='tight')

def integrate_aneuploidy_statistics(data_type):


    types = sorted(get_cancer_types_in_model(data_type))
    fieldnames = ['Type'] + get_chromosome_arms()
    threshold = get_qvalue_threshold(data_type)
    with open("../article_results/Figure 1/panel A _ cancer_arm_aneuploidy.csv", 'w',newline='') as csvfile:

        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for t in types:
            read = pd.read_csv("../Datasets/GISTIC_aneuploidy_stats_per_cancer/"+t+".csv")
            write_row = {}
            write_row['Type'] = t

            for rowNum,row in read.iterrows():
                arm = row['Arm']
                write_row[arm] = get_aneuploidy_status(row['Amp Q value'],row['Del Q value'],row['Amp Frequency'],row['Del Frequency'],threshold)

                
            writer.writerow(write_row)

    results_as_heatmap()

