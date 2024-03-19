import pandas as pd
from scipy import stats

def correlation_between_fearues():
    data_to_save = {'feature':[]}
    for col in columns:
        data_to_save[col] = []
    for i in range(0,len(columns)):
        f1 = columns[i]
        data_to_save['feature'].append(f1)
        for j in range (0,len(columns)):
            f2 = columns[j]
            data_list_f1 = list(Dataset[f1])
            data_list_f2 = list(Dataset[f2])
            for i in range(len(data_list_f1)-1,-1,-1):
                if data_list_f1[i] != data_list_f1[i] or data_list_f2[i] != data_list_f2[i]:
                    del data_list_f1[i]
                    del data_list_f2[i]

            r,p = stats.spearmanr(data_list_f1,data_list_f2)
            data_to_save[f2].append(r)
    pd.DataFrame.from_dict(data_to_save).to_csv('feature_correlation_results.csv',index = False)


Dataset = pd.read_csv("../../../article_results/Data/Model_dataset.csv")
columns = Dataset.columns.tolist()
columns.remove('Arm')
columns.remove('Type')
columns.remove('Label')
correlation_between_fearues()
