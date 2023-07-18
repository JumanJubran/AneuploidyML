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
        sample_linage[row['CCLE_Name']] = row['lineage']
        if row['lineage'] not in linage_samples:
            linage_samples[row['lineage']] = [row['CCLE_Name']]
        else:
            linage_samples[row['lineage']].append(row['CCLE_Name'])
    return linage_samples, sample_linage


def get_fieldnames(linage_samples):
    '''
    This function returns the fieldnames of the rnai features
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
    del fieldnames["fibroblast"]
    del fieldnames["thyroid"]
    del fieldnames["peripheral_nervous_system"]
    del fieldnames["bile_duct"]
    return fieldnames


def associate_names():
    '''
    This function returns a dictionary that associate each RNAi sample with DepMap sample
    :return:
    '''
    info = pd.read_csv("../Datasets/CRISPR data/internal-21q1_v30-sample-info.csv")
    DepMap_CCLE_names = {}
    for rowNum, row in info.iterrows():
        row = row.to_dict()
        DepMap_CCLE_names[row['DepMap_ID']] = row['CCLE_Name']
    return DepMap_CCLE_names


def get_gene_neutrals(gene_chromosome, neutral_samples, sample_linage, DepMap_CCLE_names):
    '''
    This function returns a dictionary that contains for each gene the RNAi scores of the neutral samples per tissue.
    :param gene_chromosome, neutral_samples, sample_linage, DepMap_CCLE_names:
    :return: gene's rnai scores of neutral samples
    '''
    RNAi = pd.read_csv("../Datasets/RNAi data/D2_Achilles_gene_dep_scores.csv")
    gene_neutrals = {}
    for rowNum, row in RNAi.iterrows():

        row = row.to_dict()
        gene = row['Unnamed: 0'].split('(')[0].strip()
        if gene in gene_chromosome:
            chromosome = gene_chromosome[gene]

            gene_neutrals[gene] = {}

            if chromosome in neutral_samples:
                for samp in neutral_samples[chromosome]:
                    CCLE_sample = DepMap_CCLE_names[samp]
                    if CCLE_sample in row:

                        linage = sample_linage[CCLE_sample]

                        if linage not in gene_neutrals[gene]:
                            gene_neutrals[gene][linage] = []
                        if row[CCLE_sample] == row[CCLE_sample]:
                            gene_neutrals[gene][linage].append(row[CCLE_sample])
    return gene_neutrals


def get_neutral_samples():
    '''
    This function returns a dictionary that contain for each chromosome arm samples that did not have CNAs (neutral samples).
    :return:
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


def add_median_rnai(Dataset, fieldnames, gene_neutrals):
    '''
    median rnai score of neutral samples per gene and tissue
    param Dataset, fieldnames, gene_neutrals:
    return: Dataset with RNAi scores
    '''
    for rowNum, row in Dataset.iterrows():
        gene = row['Gene_name']
        if gene in gene_neutrals:
            gene_rnai_scores = gene_neutrals[gene]
            for linage in fieldnames:
                if linage in gene_rnai_scores:
                    if len(gene_rnai_scores[linage]) != 0:
                        fieldnames[linage].append(stats.median(gene_rnai_scores[linage]))
                    else:
                        fieldnames[linage].append('')
                else:
                    fieldnames[linage].append('')
        else:
            for linage in fieldnames:
                fieldnames[linage].append('')

    for linage in fieldnames:

        Dataset[linage.replace('_','-') + '_(neutral) median rnai score'] = fieldnames[linage]
    return Dataset



def RNAi_Main(Dataset):
    linage_samples, sample_linage = get_linage_info()
    gene_chromosome, chromosome_genes = get_genes_chromosome_arm()  # This function from get_info file
    neutral_samples = get_neutral_samples()
    DepMap_CCLE_names = associate_names()
    gene_neutrals = get_gene_neutrals(gene_chromosome, neutral_samples, sample_linage, DepMap_CCLE_names)
    fieldnames = get_fieldnames(linage_samples)
    Dataset = add_median_rnai(Dataset, fieldnames, gene_neutrals)
    Dataset.to_csv('Processed datasets/5_TRACE_Dataset_with_RNAi.csv', index=False)
    return Dataset











