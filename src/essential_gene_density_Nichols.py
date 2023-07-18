from get_info import *


def calculate_density():

    data = pd.read_csv("../Datasets/Essential gene list from Will Gibson.txt",delimiter = '\t')
    arms = []
    gene_arms, _ = get_genes_chromosome_arm()
    for rowNum, row in data.iterrows():
        g = row['gene']

        if row['gene'].strip() in gene_arms:
            arms.append(gene_arms[row['gene']].strip())
        else:
            arms.append('')

    data['Arm'] = arms
    arms_count = {}
    arms_essentials = {}
    for rowNum, row in data.iterrows():
        arm = row['Arm']
        if arm not in arms_count:
            arms_count[arm] = 0
            arms_essentials[arm] = 0
        arms_count[arm] += 1
        if row['essentiality_probability'] >= 0.8:
            arms_essentials[arm] += 1
    i = 0

    density_dict = {}
    for arm in arms_count:
        if arm != '' and arm != '21p':
            density_dict[arm] = arms_essentials[arm]/arms_count[arm]
    return density_dict

def add_Essential_gene_density_Nichols(Dataset):

    densities = calculate_density()
    dens = []
    
    for rowNum, row in Dataset.iterrows():
        dens.append(densities[row['Arm']])
    
    Dataset['Essentiality density (Nichols)'] = dens
    return Dataset

    
   

