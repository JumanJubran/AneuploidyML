import statistics as stat

from get_info import *

def calculate_median_top(Dataset):

    chromosome_hierarchical_info = {}   ## Building set of values based on chromosome and feature  dict[chromosome][feature] = [values]
    chromosome_hierarchical_gene_info = {}  ## Building set of values based on chromosome, feature and gene  dict[chromosome][feature][genes] = value 
    chromosome_genes_count = {} ## Number of genes on chromosome
    
    for rowNum, row in Dataset.iterrows(): ## Iteration on z-score dataset
        row = row.to_dict()
        
        ### Building the dataset above (chromosome_hierarchical_info,chromosome_hierarchical_gene_info,chromosome_genes_count)
        chromosome_arm = row['Chromosome_arm']
        if chromosome_arm not in chromosome_hierarchical_info:
            chromosome_hierarchical_info[chromosome_arm] = {}
            chromosome_hierarchical_gene_info[chromosome_arm] = {}
            chromosome_genes_count[chromosome_arm] = 0
            for feature in row:
                if feature not in ['Chromosome_arm','Gene_name','Gene_ID']:
                    chromosome_hierarchical_info[chromosome_arm][feature] = []
                    chromosome_hierarchical_gene_info[chromosome_arm][feature] = {}
        
        chromosome_genes_count[chromosome_arm] += 1
        for feature in row:
            if feature not in ['Chromosome_arm','Gene_name','Gene_ID']:
                if row[feature] == row[feature] and row[feature]!='':
                    chromosome_hierarchical_info[chromosome_arm][feature].append(row[feature])
                    chromosome_hierarchical_gene_info[chromosome_arm][feature][row['Gene_name']] = row[feature]
                    

    #get_top_10_genes(chromosome_hierarchical_gene_info, z_scores, path) ##Getting top 10% genes of each feature based on z-scores

    median_top_10_per = {}
    median_top_10_per['chromosome arm'] = []
    
    ## Calculating the median value of top 10% genes 
    for chromosome in chromosome_hierarchical_info:
        if chromosome != '21p':
            median_top_10_per['chromosome arm'].append(chromosome)
            for feature in chromosome_hierarchical_info[chromosome]:
                
               
                if feature not in median_top_10_per:
                    median_top_10_per[feature] = []
                
                sorted_vals = sorted(chromosome_hierarchical_info[chromosome][feature])
                
                ten_precent = round(len(sorted_vals)*0.1)
             
                if len(sorted_vals) != 0:
                    median_top_10_per[feature].append(stat.median(sorted_vals[-ten_precent:]))
         
             
                else:
                    median_top_10_per[feature].append('')
       
    Dataset = pd.DataFrame.from_dict(median_top_10_per)
    Dataset.to_csv('Processed datasets/7_median_top_10_percent.csv', index = False)
    return Dataset
    

