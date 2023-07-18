from build_aneuploidy_table_GISTIC import *
from TRACE_preprocessing import *
from TCGA_cancer_expression import *
from RNAi_Main import *
from calculate_z_score import *
from calcuate_median_top import *
from construct_dataset import *

def dataset_processing_main():

    print('---------- Assigning aneuploidy status based on GISTIC 2.0 ----------')
    integrate_aneuploidy_statistics('TCGA')                                         #This function is from build_aneuploidy_table_GISTIC file


    print('----------TRACE Preprocessing----------')
    Dataset_genes = TRACE_preprocessing_Main()                                      #This function is from TRACE_preprocessing file


    print('----------TCGA Calculations----------')
    Dataset_genes = TCGA_Main(Dataset_genes)                                        #This function is from median_expression_per_gene file

    print('----------CRISPR Calculations----------')
    Dataset_genes = Crispr_Main(Dataset_genes)                                      #This function is from Crispr_main file

    print('----------RNAi Calculations----------')
    Dataset_genes = RNAi_Main(Dataset_genes)                                        #This function is from RNAi_main file


    print('----------Calculate z-score of genes----------')
    Dataset_genes_zscore = calculate_z_score(Dataset_genes)                        #This function is from calculate_z_score

    
    print('----------Calculate top 10% genes per feature nd tissue----------')
    Dataset_arm = calculate_median_top(Dataset_genes_zscore)                       #This function is from calcuate_median_top
    
    print('----------Construct cancer type - chromosome arm dataset for ML learning----------')
    construct_Dataset_main(Dataset_arm, Dataset_genes_zscore, Dataset_genes)

