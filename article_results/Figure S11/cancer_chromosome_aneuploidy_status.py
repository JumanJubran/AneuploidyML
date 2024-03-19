import csv
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


from get_info import *

def check_if_all_significant(arr):
    for i in arr:
        if i >= 0.05:
            return False
    return True


def whole_chromosome_labels():
    types = sorted(get_cancer_types_in_model('TCGA'))

    dict_to_write = {'Type':[]}
    for i in range(1,23):
        dict_to_write[str(i)] = []


    for cType in types:
        data = pd.read_csv('../../Datasets/GISTIC_aneuploidy_stats_per_cancer/'+cType+'.csv')
        arm_q_vals = {}
        for rowNum, row in data.iterrows():
            arm = row['Arm'].replace('p','').replace('q','')
            if arm not in arm_q_vals:
                arm_q_vals[arm] = [[],[]]
            arm_q_vals[arm][0].append(row['Amp Q value'])
            arm_q_vals[arm][1].append(row['Del Q value'])

        dict_to_write['Type'].append(cType)
        for arm in arm_q_vals:
            amp_flag = check_if_all_significant(arm_q_vals[arm][0])
            del_flag = check_if_all_significant(arm_q_vals[arm][1])

            if amp_flag == True and del_flag == False:
                dict_to_write[arm].append(1)
            if amp_flag == False and del_flag == True:
                dict_to_write[arm].append(-1)
            if amp_flag == False and del_flag == False:
                dict_to_write[arm].append(0)
            if amp_flag == True and del_flag == True:
                print(cType, arm)

    pd.DataFrame.from_dict(dict_to_write).to_csv('cancer_aneuploidy.csv', index = False)
