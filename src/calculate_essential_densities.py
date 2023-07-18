from Crispr_main import *


def calculate_densities(chromosome_genes, Dataset_genes):
    '''
    This function creats a csv file 'essential genes density.csv' that contain fro each chromosome arm the density of essential gene for each tissue
    :param chromosome_genes, Dataset_genes
    :return:
    '''

    essential_density = {'Chromosome_arm':[]}
    for item in Dataset_genes.columns.tolist():
        if 'crispr' in item:
            tissue = item.split('_')[0]
            essential_density[tissue] = []


    for chromosome in chromosome_genes:
        if chromosome != '21p':
            essential_density['Chromosome_arm'].append(chromosome)
            Dataset_genes_chrom = Dataset_genes.loc[Dataset_genes['Chromosome_arm']==chromosome]
            for tissue in essential_density:
                if tissue!='Chromosome_arm':
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
    essential_density_dataframe.to_csv('Processed datasets/essential_gene_density.csv', index = False)
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
        essential_density_dict[row['Chromosome_arm']] = row

    for rowNum, row in Dataset.iterrows():
        disease_tissue = []
        unaffected_tissue = []
        row = row.to_dict()
        cancer_type = row['Type']
        c_arm = row['Arm']

        if cancer_type in cancer_tissue_dict:
            d_tissue = cancer_tissue_dict[cancer_type]

            densities = essential_density_dict[c_arm]
            for tissue in densities:
                if tissue != 'Chromosome_arm':
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
            mean = stats.median(disease_tissue + unaffected_tissue)
            stdv = stats.stdev(disease_tissue + unaffected_tissue)
            z_score = (chromosome_score - mean) / stdv
            density_feature.append(z_score)

        else:
            density_feature.append('')

    Dataset['essential genes density'] = density_feature
    return Dataset

def crispr_densities_main(Dataset_genes, arm_cancer_dataset):

    gene_chromosome,chromosome_genes = get_genes_chromosome_arm()                                   #This function from get_info file
    essential_density = calculate_densities(chromosome_genes, Dataset_genes)
    arm_cancer_dataset = add_feature_to_dataset(arm_cancer_dataset, essential_density)
    return arm_cancer_dataset

    

    

