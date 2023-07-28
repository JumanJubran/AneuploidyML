from get_info import *

def get_fieldnames(linage_samples):
    fieldnames = {'Gene name': [], 'Arm': []}
    for linage in linage_samples:
        if linage not in ["engineered lung", "engineered prostate", "engineered ovary", "engineered breast",
                          "engineered bone", "engineered central nervous system",
                          "unknown", "adrenal cortex", "thymus", "engineered bile duct", "engineered kidney", "eye",
                          "cervix", "embryo", "epidermoid_carcinoma",
                          "engineered", "engineered blood"]:
            fieldnames[linage.replace('_',' ')] = []

    return fieldnames


def get_linage_info():
    Crispr_samples = pd.read_csv("../Datasets/CRISPR data/internal-21q1_v30-sample-info.csv")

    linage_samples = {}
    sample_linage = {}

    for rowNum, row in Crispr_samples.iterrows():
        linage = row['lineage'].replace('_',' ')
        sample_linage[row['DepMap_ID']] = linage
        if linage not in linage_samples:
            linage_samples[linage] = [row['DepMap_ID']]
        else:
            linage_samples[linage].append(row['DepMap_ID'])
    return linage_samples, sample_linage


def get_crisprs_per_gene_linage(crispr, gene_chromosome, sample_linage):
    gene_crispr = {}
    for rowNum, row in crispr.iterrows():

        row = row.to_dict()
        gene = row['Unnamed: 0']
        if gene in gene_chromosome:
            chromosome = gene_chromosome[gene]

            gene_crispr[gene] = {}

            for samp in row:
                if samp != 'Unnamed: 0':
                    linage = sample_linage[samp]

                    if linage not in gene_crispr[gene]:
                        gene_crispr[gene][linage] = []
                    gene_crispr[gene][linage].append(row[samp])
    return gene_crispr

def calculate_min_crispr():
    crispr = pd.read_csv("../Datasets/CRISPR data/internal-21q1_v30-achilles-gene-effect-T.csv")

    gene_chromosome, chromosome_genes = get_genes_chromosome_arm()  # This function from get_info file

    linage_samples, sample_linage = get_linage_info()
    gene_crispr = get_crisprs_per_gene_linage(crispr, gene_chromosome, sample_linage)
    gene_chromosome, chromosome_genes = get_genes_chromosome_arm()

    dict_to_write = get_fieldnames(linage_samples)

    for gene in gene_crispr:

        dict_to_write['Gene name'].append(gene)

        if gene in gene_chromosome:
            dict_to_write['Arm'].append(gene_chromosome[gene])
        else:
            dict_to_write['Arm'].append('')

        for linage in dict_to_write:
            if linage not in ['Gene name', 'Arm']:
                if linage in gene_crispr[gene]:
                    dict_to_write[linage].append(min(gene_crispr[gene][linage]))
                else:
                    dict_to_write[linage].append('')

    pd.DataFrame.from_dict(dict_to_write).to_csv('Processed datasets/crispr_min_score_all_samples.csv',index=False)