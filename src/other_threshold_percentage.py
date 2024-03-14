from construct_dataset import *
from calcuate_median_top import *
import pandas as pd

for threshold in [0.01,0.05,0.15,0.2]:
    Dataset_gene_zscore = pd.read_csv('Processed datasets/6_Genes z score.csv')
    Dataset_genes = pd.read_csv('Processed datasets/5_TRACE_Dataset_with_RNAi.csv')

    percent = str(int(threshold*100))
    Dataset_arm = calculate_median_top(Dataset_gene_zscore,threshold)                       #This function is from calcuate_median_top
    Dataset_arm.to_csv('../article_results/Figure S7/'+percent+' percent/median_top_'+percent+'_percent_new.csv', index=False)


    arm_cancer_dataset = construct_Dataset_main(Dataset_arm, Dataset_gene_zscore, Dataset_genes, 'TCGA')

    arm_cancer_dataset.to_csv('../article_results/Figure S7/'+percent+' percent/Arm_Cancer_Dataset_top_'+percent+'_per_new.csv', index=False)
