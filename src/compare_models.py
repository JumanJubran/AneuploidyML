import pandas as pd
import statistics
import numpy as np
import csv


#Models:
from sklearn.multiclass import OneVsRestClassifier
import xgboost as xgb
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import BaggingClassifier
from sklearn.neural_network import MLPClassifier

#TenFold:
from sklearn.model_selection import StratifiedKFold


#Plotting:
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import auc


'''
Tested 5 ML models: XGboost, Gradient Boosting, Random Forest, Bagging and Logistic regression
'''
def get_models():

    return {
    'XGBoost': xgb.XGBClassifier( eval_metric='logloss'),
    'Random Forest': RandomForestClassifier(),
    'Logistic Regression': LogisticRegression(solver='liblinear'),
    'Gradient Boosting': GradientBoostingClassifier(),
    'Bagging': BaggingClassifier(),
    }


'''
pointing colors per each model for the curves plots
'''
def get_models_color():

    return {
    'XGBoost': '#66c2a5',
    'Random Forest': '#fc8d62',
    'Logistic Regression': '#8da0cb',
    'Gradient Boosting': '#e78ac3',
    'Bagging': '#a6d854',
    }


'''
Computing false positive rate (fpr), true positive rate (tpr), recall and precision
'''
def compute_plot_vectors(tests, predictions):
    
   
    fpr, tpr, _ = roc_curve(tests, predictions)
    roc_auc = auc(fpr, tpr)
   
    prec, recall, _ = precision_recall_curve(tests, predictions)
    PR_auc = auc(recall,prec)

    return fpr,tpr,roc_auc,prec,recall,PR_auc


'''
plotting ROC and PR curves of all 5 models
'''
def plot_curves(model_roc_scores,roc_median,model_PR_scores,pr_median,file_name, threshold):
    '''
        plot Roc Curve and percision-recall curve.
    '''
    models_dict = get_models()
    models_colors = get_models_color()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 8))
    for m in models_dict:
        model_roc_median = roc_median[m]
        model_pr_median = pr_median[m]
        
        lw=3
        
        ax1.plot(model_roc_scores[m][model_roc_median][0], model_roc_scores[m][model_roc_median][1], lw=lw, color = models_colors[m], label = m )
        ax2.plot(model_PR_scores[m][model_pr_median][0], model_PR_scores[m][model_pr_median][1], lw=lw, color = models_colors[m], label=m )


        ax1.plot([0, 1], [0, 1], color='black', lw=lw, linestyle='--')
        ax1.set_xlim([0.0, 1.0])
        ax1.set_ylim([0.0, 1.05])
        ax1.set_xlabel('False Positive Rate')
        ax1.set_ylabel('True Positive Rate')
        ax1.set_title('ROC curve')
        ax1.legend(loc="lower right")
     
        ax2.plot([0, 1], [threshold, threshold], color='black', lw=lw, linestyle='--')   
        ax2.set_xlim([0.0, 1.0])
        ax2.set_ylim([0.0, 1.05])
        ax2.set_xlabel('Recall')
        ax2.set_ylabel('Percision')
        ax2.set_title('Percision-Recall curve')
        ax2.legend(loc="lower right")
    plt.savefig('../article_results/Figure S3/'+ file_name + "__ROC_PR_Curves.png",bbox_inches='tight')



'''
Finding the median fold
'''
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


'''
Getting the information of the model's dataset
'''
def get_data_info(c1,c2,counter):
    total_count = counter['total']
    if c1 == 'Loss':
        class_to_drop = 1
        labels = {-1:1, 0:0, 1:0}
        threshold = counter[-1]/counter['total']
        if c2 == 'Neutral':
            
            total_count = counter[-1]+counter[0]
            threshold = counter[-1]/total_count
    else:
        class_to_drop = -1
        labels = {-1:0, 0:0, 1:1}
        threshold = counter[1]/counter['total']
        if c2 == 'Neutral':
            total_count = counter[1]+counter[0]
            threshold = counter[1]/total_count
            


  
                
    return labels, class_to_drop, threshold, total_count


'''
process the data based on the model. there are 4 possible models, including: 
Gain vs. Neutral, 
Loss vs. Neutral, 
Gain vs. Rest and 
Loss vs. Rest.
'''
def process_dataset_per_model(data,c1,c2):
    num_observation_per_class = {x: list(data['Label']).count(x) for x in list(data['Label'])}
    num_observation_per_class['total'] = len(list(data['Label']))

    labels, class_to_drop, threshold, total_count = get_data_info(c1, c2, num_observation_per_class)

    if c2 == 'Neutral':
        for rowNum, row in data.iterrows():
            if row['Label'] == class_to_drop:
                data = data.drop(rowNum)

    data = data.replace({'Label': labels})
    data = data.set_axis(range(0, total_count))

    return data, threshold


'''
Fill null values with the median value of the corresponding feature
'''
def null_to_median(Dataset):
    feature_median = {}

    for column in Dataset.columns.tolist():
        values_median = []
        values = Dataset[column].values.tolist()
        for val in values:
            if val == val:
                values_median.append(val)
        feature_median[column] = statistics.median(values_median)
    for rowNum, row in Dataset.iterrows():
        row = row.to_dict()
        for item in row:
            if row[item] != row[item]:
                Dataset.at[rowNum, item] = feature_median[item]
    return Dataset


'''
Save the results as csv output file
'''
def get_results_as_CSV_file(X, y, models_dict, file_name, threshold,model_roc_scores, model_PR_scores):
    flag = 0
    results_output_file = {'Fold':['1','2','3','4','5','6','7','8','9','10','Average']}
    for m in model_roc_scores:
        roc_vals = model_roc_scores[m].keys()
        results_output_file[m + ' ROC'] = list(roc_vals)
        results_output_file[m + ' ROC'].append(statistics.mean(list(roc_vals)))

        pr_vals = model_PR_scores[m].keys()
        results_output_file[m + ' PR'] = list(pr_vals)
        results_output_file[m + ' PR'].append(statistics.mean(list(pr_vals)))
    for i in results_output_file:
        if len(results_output_file[i]) < 11:
            flag = 1
            results_output_file[i].append('')
            #methods_comparision_by_ten_fold_validation(X, y, models_dict, file_name, threshold)

    if flag == 0:
        pd.DataFrame.from_dict(results_output_file).to_csv('../article_results/Figure S3/'+file_name + '__Ten_fold_results.csv', index=False)


'''
calculating the median roc and pr, per method
'''
def find_roc_pr_median(model_roc_scores,model_PR_scores):
    roc_median = {}
    pr_median = {}
    for m in model_roc_scores:
        roc_median[m] = {}
        pr_median[m] = {}

        roc_vals = model_roc_scores[m].keys()
        roc_med = statistics.mean(roc_vals)
        roc_median[m] = find_nearest(list(roc_vals), roc_med)

        pr_vals = model_PR_scores[m].keys()
        pr_med = statistics.mean(pr_vals)
        pr_median[m] = find_nearest(list(pr_vals), pr_med)
    return roc_median, pr_median



'''
comparing between different ML methods
'''
def methods_comparision_by_ten_fold_validation(X, y, models_dict, file_name, threshold):

    model_roc_scores = {}
    model_PR_scores = {}

    ### Ten fold cross validation by StratifiedKFold
    kf = StratifiedKFold(n_splits=10, shuffle=True, random_state=1001)

    ### Each fold are divided into train and test
    for train, test in kf.split(X, y):

        train_index = train.tolist()
        test_index = test.tolist()

        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y.iloc[train_index], y.iloc[test_index]

        ### Test all ML methods
        for model_name in models_dict:

            model = models_dict[model_name]

            predictions = model.fit(X_train, y_train).predict_proba(X_test)[:, 1]

            # Compute ROC curve and ROC area for each class
            fpr, tpr, roc_auc, prec, recall, PR_auc = compute_plot_vectors(np.asarray(y_test.tolist()), predictions)

            if model_name not in model_roc_scores:
                model_roc_scores[model_name] = {}
                model_PR_scores[model_name] = {}

            model_roc_scores[model_name][roc_auc] = [fpr, tpr]
            model_PR_scores[model_name][PR_auc] = [recall, prec]

    ### Save results
    get_results_as_CSV_file(X, y, models_dict, file_name, threshold,model_roc_scores, model_PR_scores)

    ### Find median fold per model
    roc_median, pr_median = find_roc_pr_median(model_roc_scores, model_PR_scores)

    # plot ROC and PR curves
    plot_curves(model_roc_scores,roc_median,model_PR_scores,pr_median,file_name, threshold)


'''
Assessing the performance of 5 different ML model per all 4 different models that includes:
Gain vs. Neutral, 
Loss vs. Neutral, 
Gain vs. Rest and 
Loss vs. Rest
'''
def compare_ML_methods_performance():
    Dataset = pd.read_csv("../article_results/Data/Model_dataset.csv")

    models_dict = get_models()


    first_class = ['Gain','Loss']
    second_class = ['Neutral','Rest']

    for c1 in first_class:
        for c2 in second_class:


            file_name = c1+' versus '+c2
            dataset_per_model, threshold = process_dataset_per_model(Dataset,c1,c2)

            dataset_per_model = dataset_per_model.drop(columns= ['Type','Arm'])

            dataset_per_model = null_to_median(dataset_per_model)

            X = dataset_per_model .loc[:, dataset_per_model .columns != 'Label']
            y = dataset_per_model ['Label']


            methods_comparision_by_ten_fold_validation(X, y, models_dict, file_name, threshold)



