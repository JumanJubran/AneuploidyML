from get_info import *
import pandas as pd
import statistics as stats


def calculate_densities(chromosome_genes, Dataset_genes):
    '''
    This function creats a csv file 'essential genes density.csv' that contain fro each chromosome arm the density of essential gene for each tissue
    :param chromosome_genes, Dataset_genes
    :return:
    '''

    essential_density = {'Chromosome':[]}
    for item in Dataset_genes.columns.tolist():
        if 'crispr' in item:
            tissue = item.split('_')[0]
            essential_density[tissue] = []


    for chromosome in chromosome_genes:

        essential_density['Chromosome'].append(chromosome)
        Dataset_genes_chrom = Dataset_genes.loc[Dataset_genes['Chromosome']==int(chromosome)]

        for tissue in essential_density:
            if tissue!='Chromosome':
                crispr_tissue = list(Dataset_genes_chrom[tissue+'_(neutral) median crispr score'])
                total_genes = 0
                essential_genes = 0
                for value in crispr_tissue:
                    if value==value and value!='':
                        if value <= -0.5:
                            essential_genes += 1
                        total_genes += 1
                if total_genes!=0:
                    essential_density[tissue].append(essential_genes/total_genes)
                else:
                    essential_density[tissue].append('')
    essential_density_dataframe = pd.DataFrame.from_dict(essential_density)
    essential_density_dataframe.to_csv('processed_datasets/essential_gene_density.csv', index = False)
    return essential_density_dataframe

def add_feature_to_dataset(Dataset, essential_density):
    '''
    This function adds the tissue specific densities to the model's features as a z-score per tissue.
    :param Dataset:
    :param ts_gene_densities:
    :return:
    '''
    essential_density_dict = {}
    cancer_tissue_dict = get_cancers_associated_tissue('TCGA', 'DepMap')  # This function from get_info file
    density_feature = []
    for rowNum, row in essential_density.iterrows():
        row = row.to_dict()
        essential_density_dict[row['Chromosome']] = row

    for rowNum, row in Dataset.iterrows():
        disease_tissue = []
        unaffected_tissue = []
        row = row.to_dict()
        cancer_type = row['Type']
        c_arm = row['Chromosome']

        if cancer_type in cancer_tissue_dict:
            d_tissue = cancer_tissue_dict[cancer_type]

            densities = essential_density_dict[str(int(c_arm))]
            for tissue in densities:
                if tissue != 'Chromosome':
                    tissue_ = tissue.replace('-', ' ').strip().lower()
                    if d_tissue in tissue_ and d_tissue != '':
                        if densities[tissue] == densities[tissue] and densities[tissue]!='':
                            disease_tissue.append(densities[tissue])
                    else:
                        if densities[tissue] == densities[tissue] and densities[tissue] != '':
                            unaffected_tissue.append(densities[tissue])

        if len(disease_tissue) != 0:

            '''
            Z-score calculation
            '''

            chromosome_score = stats.median(disease_tissue)
            mean = stats.mean(disease_tissue + unaffected_tissue)
            stdv = stats.stdev(disease_tissue + unaffected_tissue)
            z_score = (chromosome_score - mean) / stdv
            density_feature.append(z_score)

        else:
            density_feature.append('')

    Dataset['essential genes density'] = density_feature
    return Dataset

def crispr_densities_main(gene_dataset,chromosome_cancer_dataset):

    gene_chromosome,chromosome_genes = get_genes_chromosome_arm()                                   #This function from get_info file

    whole_chromosome_genes = {}
    for chromosome_arm in chromosome_genes:
        chromosome = chromosome_arm.replace('q','').replace('p','')
        if chromosome not in whole_chromosome_genes:
            whole_chromosome_genes[chromosome] = []
        whole_chromosome_genes[chromosome] += chromosome_genes[chromosome_arm]




    essential_density = calculate_densities(whole_chromosome_genes, gene_dataset)

    chromosome_cancer_dataset = add_feature_to_dataset(chromosome_cancer_dataset, essential_density)
    
    return chromosome_cancer_dataset


    

