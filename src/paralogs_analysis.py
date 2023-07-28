from crispr_cal_min import *
from get_info import *
import scipy.stats as stats
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
def get_aneuploidy_label():
    data = pd.read_csv("../article_results/Figure 1/panel A _ cancer_arm_aneuploidy.csv")
    aneuploidy_labels = {}
    for rowNum, row in data.iterrows():
        row = row.to_dict()
        for item in row:
            if item != 'Type':
                aneuploidy_labels[(row['Type'],item)] = row[item]
    return aneuploidy_labels
def paralogous_pairs():
    paralog_data = pd.read_csv( "../Datasets/median-highest identity paralog.csv")
    paralogs = {}
    for rowNum, row in paralog_data.iterrows():
        paralogs[row['gene Description']] = row['paralog Description']
    return paralogs

def get_essential_scores(gene, tissue,crispr):
    cris = crispr.loc[crispr['Gene name'] == gene].to_dict('records')
    if len(cris) > 0:
        return cris[0][tissue]
    else:
        return ''

def add_essen_score(essen_gained_paralog,non_essen_gained_paralog,intermediate_gained_paralog,essen_score):
    essential_threshold = -0.5
    non_essential_threshold = -0.3
    if essen_score <= essential_threshold:
        essen_gained_paralog.append(essen_score)
    elif essen_score >= non_essential_threshold:
        non_essen_gained_paralog.append(essen_score)
    else:
        intermediate_gained_paralog.append(essen_score)
    return essen_gained_paralog,non_essen_gained_paralog,intermediate_gained_paralog

def histogram(x, y, color1, color2, label1, label2, name):
    fig = plt.figure()

    bins = np.linspace(-3, 0.5, 50)

    plt.hist([x, y], bins, label=[label1, label2], weights=[np.ones(len(x)) / len(x), np.ones(len(y)) / len(y)],
             color=[color1, color2])

    # plt.hist(x, bins=100, alpha=0.5, label=label1,weights = np.ones(len(x)) / len(x), color = color1, histtype="step")
    # plt.hist(y, bins=100, alpha=0.5, label=label2,weights = np.ones(len(y)) / len(y), color = color2, histtype="step")

    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.legend(loc='upper left')
    plt.savefig('../article_results/Figure S8/paralogs histogram.png')

def lost_paralogs_analysis():
    gain_color = '#BF0000'
    loss_color = '#4472C4'

    crispr = pd.read_csv('Processed datasets/crispr_min_score_all_samples.csv')

    aneuploidy_labels = get_aneuploidy_label()  ### get aneuploidy label per pair of arm and cancer
    cancer_type = get_cancer_types_in_model('TCGA')  ### get cancers in the analysis
    chromosome_arms = get_chromosome_arms()  ### get chromosome arms
    gene_chromosome, chromosome_genes = get_genes_chromosome_arm()  ### get pairs of gene and arm
    paralogs = paralogous_pairs()  ### get paralogous pairs
    cancer_tissues = get_cancers_associated_tissue('TCGA', 'DepMap')  ### get per cancer its tissue of origin

    essen_gained_paralog, intermediate_gained_paralog, non_essen_gained_paralog  = [],[],[]
    essen_lost_paralog, intermediate_lost_paralog,non_essen_lost_paralog  = [],[],[]
    essen_neutral_paralog, intermediate_neutral_paralog, non_essen_neutral_paralog  = [],[],[]
    gained_essen = []
    lost_essen = []

    for cancer in cancer_type:  ### for each cancer

        if cancer in cancer_tissues:

            for arm in chromosome_arms:  ### for each arm

                if aneuploidy_labels[(cancer, arm)] == -1:  ### lost genes only
                    tissue = cancer_tissues[cancer]

                    for gene in chromosome_genes[arm]:

                        if gene in paralogs:
                            paralog = paralogs[gene]
                            if paralog in gene_chromosome:
                                para_arm = gene_chromosome[paralog]

                                essen_score = get_essential_scores(gene, tissue,crispr)
                                if essen_score == essen_score and essen_score != '':

                                    #### Gained paralogs
                                    if aneuploidy_labels[(cancer, para_arm)] == 1:
                                        gained_essen.append(essen_score)
                                        essen_gained_paralog, non_essen_gained_paralog, intermediate_gained_paralog = add_essen_score(
                                            essen_gained_paralog, non_essen_gained_paralog, intermediate_gained_paralog,
                                            essen_score)
                                    #### Lost paralogs
                                    elif aneuploidy_labels[(cancer, para_arm)] == -1:
                                        lost_essen.append(essen_score)
                                        essen_lost_paralog, non_essen_lost_paralog, intermediate_lost_paralog = add_essen_score(
                                            essen_lost_paralog, non_essen_lost_paralog, intermediate_lost_paralog,
                                            essen_score)
                                    #### Neutral paralogs
                                    else:
                                        essen_neutral_paralog, non_essen_neutral_paralog, intermediate_neutral_paralog = add_essen_score(
                                            essen_neutral_paralog, non_essen_neutral_paralog,
                                            intermediate_neutral_paralog, essen_score)

    results_dict = {'Paralog type': ['Gained paralog', 'Neutral paralog', 'Lost paralog'],
                    'Essential': [len(essen_gained_paralog), len(essen_neutral_paralog), len(essen_lost_paralog)],
                    'Intermediate': [len(intermediate_gained_paralog), len(intermediate_neutral_paralog),
                                     len(intermediate_lost_paralog)],
                    'Non-essential': [len(non_essen_gained_paralog), len(non_essen_neutral_paralog),
                                      len(non_essen_lost_paralog)]}

    pd.DataFrame.from_dict(results_dict).to_csv('../article_results/Figure S8/Panel B.csv')

    histogram(gained_essen,lost_essen,gain_color,loss_color,'gained paralogs','lost_paralogs',"gain vs lost of lost genes")
    print(stats.kstest(lost_essen, gained_essen, alternative='less'))


def analyze_lost_paralogs():
    calculate_min_crispr()
    lost_paralogs_analysis()
