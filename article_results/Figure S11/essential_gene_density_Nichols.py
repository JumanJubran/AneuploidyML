from get_info import *


def calculate_density():

    data = pd.read_csv("../../Datasets/Essential gene list from Will Gibson.txt",delimiter = '\t')
    chromosome = []
    gene_arms, _ = get_genes_chromosome_arm()
    for g in gene_arms:
        gene_arms[g] = gene_arms[g].replace('q','').replace('p','')

    for rowNum, row in data.iterrows():
        g = row['gene']

        if row['gene'].strip() in gene_arms:
            chromosome.append(gene_arms[row['gene']].strip())
        else:
            chromosome.append('')

    data['Chromosome'] = chromosome
    chromosomes_count = {}
    chromosomes_essentials = {}
    for rowNum, row in data.iterrows():
        chromosome = row['Chromosome']
        if chromosome not in chromosomes_count:
            chromosomes_count[chromosome] = 0
            chromosomes_essentials[chromosome] = 0
        chromosomes_count[chromosome] += 1
        if row['essentiality_probability'] >= 0.8:
            chromosomes_essentials[chromosome] += 1
    i = 0

    density_dict = {}
    for chromosome in chromosomes_count:

        density_dict[chromosome] = chromosomes_essentials[chromosome]/chromosomes_count[chromosome]
    return density_dict

def add_Essential_gene_density_Nichols(Dataset):

    densities = calculate_density()
    dens = []
    
    for rowNum, row in Dataset.iterrows():
        dens.append(densities[str(int(row['Chromosome']))])
    
    Dataset['Essentiality density (Nichols)'] = dens
    return Dataset

    
   

