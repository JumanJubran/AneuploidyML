from calculate_z_score import *
from calcuate_median_top import *
from build_aneuploidy_table_GISTIC import *
from construct_dataset import *
from get_info import *
from CCL_model_analysis import *


def get_tissue(sample):
    sample = sample.split('_')
    t = ''
    for i in range(1,len(sample)):
        t += sample[i] + ' '
    t = t.strip().lower()
    
    return t

def add_fieldnames(data_dict,samples):

    for rowNum, row in samples.iterrows():
        row = row.to_dict()
        tissue = get_tissue(row['CCLE_ID'])
        tissue_feature = tissue+'_(neutral)_CCL_expression'
        if tissue_feature not in data_dict:
            data_dict[tissue_feature] = []
    return data_dict
        
def neutral_gene_median_expression():
    samples = pd.read_csv("../Datasets/CCL data/Cohen-Sharir_Table_S1.csv")
    CNA_chromosomes = {}
    depMap_to_CCLE = {}
    #finding neutral chromosome arms per sample
    CCLE_tissue = {}
    for rowNum, row in samples.iterrows():
        
        row = row.to_dict()
        CNA_chromosomes [row['CCLE_ID']] = []
        depMap_to_CCLE[row['CCLE_ID']] = row['DepMap_ID']
        CCLE_tissue[row['CCLE_ID']] = get_tissue(row['CCLE_ID'])
        for item in row:
            
            if item not in ['DepMap_ID','CCLE_ID','ploidy']:
                if row[item] == 0:
                    CNA_chromosomes[row['CCLE_ID']].append(item)

    #gene_chromosome_dict holds for each gene its location on which chromosome arm.

    gene_chromosome_dict, chromosome_genes = get_genes_id_chromosome_arm()
    gene_id_to_name, gene_name_to_id = get_genes_name()
    expression_dict = {'Gene_ID':[], 'Gene_name':[], 'Chromosome_arm':[]}
    
    
    
    CCL_expression = pd.read_csv("../Datasets/CCL data/CCLE_RNAseq_genes_rpkm_20180929.csv")
    expression_dict = add_fieldnames(expression_dict,samples)

    for rowNum, row in CCL_expression.iterrows():

        row = row.to_dict()
        expressions_by_tissue = {}
        gene_id = row['Name'][:15]

        if gene_id in gene_chromosome_dict:
            chrom = gene_chromosome_dict[gene_id]

            for sample in row:

                if sample not in ['Name','Description']:
                    if sample in depMap_to_CCLE:

                        tissue = CCLE_tissue[sample].replace('_','-')                                       #Get the tissue of the sample

                        if sample in CNA_chromosomes:                                       #Check if the sample is the aneuploidy table
                            if chrom in CNA_chromosomes[sample]:                            #Check if chromosome is neutral of this specific sample
                                if tissue not in expressions_by_tissue:                     #If tissue is not in the dictionary add it
                                    expressions_by_tissue[tissue] = []
                                expressions_by_tissue[tissue].append(row[sample])  #if the chromosome is neutral for this specific sample add the expression of the gene in the sample.


        if len(expressions_by_tissue)!=0:
            expression_dict['Gene_ID'].append(gene_id)
            expression_dict['Chromosome_arm'].append(chrom)
            expression_dict['Gene_name'].append(gene_id_to_name[gene_id])
            for item in expression_dict:
                if item not in ['Gene_ID','Chromosome_arm','Gene_name']:
                    item_s = item.split('_')[0]
                    if item_s in expressions_by_tissue:
                        expression_dict[item].append(statistics.median(expressions_by_tissue[item_s]))
                    else:
                        expression_dict[item].append('')


    expression_df = pd.DataFrame.from_dict(expression_dict)
    return expression_df


def add_ccl_expression_to_dataset(arm_dataset,expression_top_median):
    features_to_add = {}
    arms = list(arm_dataset['chromosome arm'])
    for arm in arms:
        arm_expressions = expression_top_median.loc[expression_top_median['chromosome arm'] == arm].to_dict('records')[0]

        for item in arm_expressions:
            if item!='chromosome arm':
                if item not in features_to_add:
                    features_to_add[item] = []
                features_to_add[item].append(arm_expressions[item])

    for item in features_to_add:
        arm_dataset[item] = features_to_add[item]
    return arm_dataset


def CCL_calculate():

    integrate_aneuploidy_statistics('CCL')
    expression_dataset = neutral_gene_median_expression()

    expression_z_scores = calculate_z_score(expression_dataset)
    expression_top_median = calculate_median_top(expression_z_scores)

    arm_dataset = pd.read_csv('Processed datasets/7_median_top_10_percent.csv')
    arm_dataset = add_ccl_expression_to_dataset(arm_dataset,expression_top_median)

    arm_dataset = construct_Dataset_main(arm_dataset, '', '', 'CCL')
    arm_dataset = arm_dataset.drop(columns = ['Development variations','Fetus development','Childhood development','Young development','TCGA expression'])
    arm_dataset.to_csv('Processed datasets/9_CCL_arm_aneuploidy.csv', index = False)

    calculate_CCL_models()



