import statistics as stats

from get_info import *


def get_linage_info():
    '''
    This function returns two dictionaries, the first contains for each linage (tissue) the samples that were taken fron it,
    the second contains for each sample the linage (tissue) that were taken from.
    :return: linage_samples, sample_linage
    '''
    Crispr_samples = pd.read_csv("../Datasets/CRISPR data/internal-21q1_v30-sample-info.csv")

    linage_samples = {}
    sample_linage = {}

    for rowNum, row in Crispr_samples.iterrows():
        sample_linage[row['DepMap_ID']] = row['lineage']

        if row['lineage'] not in linage_samples:
            linage_samples[row['lineage']] = [row['DepMap_ID']]
        else:
            linage_samples[row['lineage']].append(row['DepMap_ID'])
    return linage_samples, sample_linage

def get_fieldnames(linage_samples):
    '''
    his function returns the fieldnames of the density file.
    :param linage_samples:
    :return:
    '''
    fieldnames = {}
    for linage in linage_samples:
        fieldnames[linage] = []
    del fieldnames["engineered_lung"]
    del fieldnames["engineered_prostate"]
    del fieldnames["engineered_ovary"]
    del fieldnames["engineered_breast"]
    del fieldnames["engineered_bone"]
    del fieldnames["engineered_central_nervous_system"]
    del fieldnames["unknown"]
    del fieldnames["adrenal_cortex"]
    del fieldnames["thymus"]
    del fieldnames["engineered_bile_duct"]
    del fieldnames["engineered_kidney"]
    del fieldnames["eye"]
    del fieldnames["cervix"]
    del fieldnames["embryo"]
    del fieldnames["epidermoid_carcinoma"]
    del fieldnames["engineered"]
    del fieldnames["engineered_blood"]
    return fieldnames

def get_gene_neutrals(gene_chromosome, neutral_samples, sample_linage):
    '''
    This function returns a dictionary that contains for each gene the crispr scores of the neutral samples per tissue.
    :param gene_chromosome, neutral_samples, sample_linage:
    :return: gene crispr score of neutral samples
    '''
    crispr = pd.read_csv("../Datasets/CRISPR data/internal-21q1_v30-achilles-gene-effect-T.csv")
    gene_neutrals = {}
    for rowNum, row in crispr.iterrows():
        row = row.to_dict()
        gene = row['Unnamed: 0']
        if gene in gene_chromosome:
            chromosome = gene_chromosome[gene]

            gene_neutrals[gene] = {}

            if chromosome in neutral_samples:
                for sample in neutral_samples[chromosome]:
                    if sample in row:
                        linage = sample_linage[sample]

                        if linage not in gene_neutrals[gene]:
                            gene_neutrals[gene][linage] = []
                        gene_neutrals[gene][linage].append(row[sample])
    return gene_neutrals

def get_neutral_samples():
    '''
    This function returns a dictionary that contain for each chromosome arm CRISPR samples that did not have CNAs (neutral samples).
    :return: chrom_neutrals
    '''
    CNAs = pd.read_csv("../Datasets/CRISPR data/Arm-level_CNAs.csv")
    chrom_neutrals = {}
    for rowNum, row in CNAs.iterrows():
        row = row.to_dict()
        for item in row:
            if item != 'Unnamed: 0':
                if item not in chrom_neutrals:
                    chrom_neutrals[item] = []
                if row[item] == 0:
                    chrom_neutrals[item].append(row['Unnamed: 0'])
    return chrom_neutrals



'''
    Function: median
    Return: NONE!
    Input files: NONE!
    Input data: fieldnames, gene_neutrals
    For each gene we calculated the median crispr score per tissue.
'''


def add_median_crispr(Dataset, fieldnames, gene_neutrals):
    '''
    median crispr score of neutral samples per gene and tissue
    param Dataset, fieldnames, gene_neutrals:
    return: Dataset with CRISPR scores
    '''
    for rowNum, row in Dataset.iterrows():
        gene = row['Gene_name']
        if gene in gene_neutrals:
            gene_crispr_scores = gene_neutrals[gene]
            for linage in fieldnames:
                if linage in gene_crispr_scores:
                    fieldnames[linage].append(stats.median(gene_crispr_scores[linage]))
                else:
                    fieldnames[linage].append('')
        else:
            for linage in fieldnames:
                fieldnames[linage].append('')

    for linage in fieldnames:
        Dataset[linage.replace('_','-') + '_(neutral) median crispr score'] = fieldnames[linage]
    return Dataset



def Crispr_Main(Dataset):
    linage_samples, sample_linage = get_linage_info()
    gene_chromosome, chromosome_genes = get_genes_chromosome_arm()  # This function from get_info file
    neutral_samples = get_neutral_samples()

    gene_neutrals = get_gene_neutrals(gene_chromosome, neutral_samples, sample_linage)
    fieldnames = get_fieldnames(linage_samples)
    Dataset = add_median_crispr(Dataset, fieldnames, gene_neutrals)
    Dataset.to_csv('Processed datasets/4_TRACE_Dataset_with_CRISPR.csv', index=False)
    return Dataset


