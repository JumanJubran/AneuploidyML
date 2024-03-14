import pandas as pd
import statistics


# Models:
import xgboost as xgb
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import BaggingClassifier

from sklearn.metrics import make_scorer, precision_score


# TenFold:
from sklearn.model_selection import StratifiedKFold


# hyperparameter Tuning
from sklearn.model_selection import RandomizedSearchCV

'''
Model's hyperparameters
'''
def get_models_param(model):
    models_params = {'XGBoost': {'eta': [0.0001, 0.001, 0.01, 0.1, 0.2, 0.3],
                                 'min_child_weight': [1, 3, 5, 10],
                                 'gamma': [0, 0.5, 1, 1.5, 2, 2.5, 3],
                                 'subsample': [0.6, 0.8, 1.0],
                                 'colsample_bytree': [0.6, 0.8, 1.0],
                                 'max_depth': [3, 4, 5, 6, 7, 10, 15],
                                 'alpha': [1, 2, 3, 4],
                                 'use_label_encoder': [False]},

                     'Random Forest': {'bootstrap': [True, False],
                                       'max_depth': [10, 20, 30, 40, None],
                                       'min_samples_leaf': [1, 2, 3, 4, 5, 10],
                                       'min_samples_split': [2, 5, 10],
                                       'n_estimators': [100, 150, 200, 300, 400]},

                     'Logistic Regression': {'penalty': ['l1', 'l2', 'none'],
                                             'C': [1, 2, 3],
                                             'solver': ['liblinear']},

                     'Gradient Boosting': {'learning_rate': [0.05, 0.1, 0.2, 0.3],
                                           'n_estimators': [100, 200, 300, 500],
                                           'subsample': [0.1, 0.3, 0.5, 0.8, 1.0],
                                           'min_samples_split': [2, 5, 10],
                                           'min_samples_leaf': [1, 3, 5, 10]},

                     'Bagging': {'n_estimators': [10, 20, 30, 50, 80, 100],
                                 'max_samples': [0.1, 0.3, 0.5, 0.8, 1.0],
                                 'max_features': [0.1, 0.3, 0.5, 0.8, 1.0],
                                 'bootstrap': [True, False]}}
    return models_params[model]

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
Getting the information of the model's dataset
'''
def get_data_info(c1, c2, counter):
    total_count = counter['total']
    if c1 == 'Loss':
        class_to_drop = 1
        labels = {-1: 1, 0: 0, 1: 0}
        threshold = counter[-1] / counter['total']
        if c2 == 'Neutral':
            total_count = counter[-1] + counter[0]
            threshold = counter[-1] / total_count
    else:
        class_to_drop = -1
        labels = {-1: 0, 0: 0, 1: 1}
        threshold = counter[1] / counter['total']
        if c2 == 'Neutral':
            total_count = counter[1] + counter[0]
            threshold = counter[1] / total_count

    return labels, class_to_drop, threshold, total_count


'''
process the data based on the model. there are 4 possible models, including: 
Gain vs. Neutral, 
Loss vs. Neutral, 
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

def add_keys(dictionary, names):
    for n in names:
        dictionary[n] = []
    return dictionary

def finding_hyperparameters(X, y,c1,c2,models_dict):
    param_comb = 200
    best_params = {}
    best_params = add_keys(best_params, models_dict)
    scorers = make_scorer(precision_score)

    for m in models_dict:
        model = models_dict[m]
        params = get_models_param(m)

        skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=1001)

        random_search = RandomizedSearchCV(model, param_distributions=params, n_iter=param_comb,
                                           scoring='average_precision', n_jobs=None, cv=skf.split(X, y), verbose=0,
                                           random_state=1001)

        random_search.fit(X, y)
        results = pd.DataFrame(random_search.cv_results_)
        print(random_search.best_score_)

        results.to_csv("../article_results/Figure S8/hyperparameter "+m+" _ "+c1+" versus "+c2+".csv", index=False)

        best_params[m].append(random_search.best_params_)
    return best_params

def parameter_tunning():

    Dataset = pd.read_csv("../article_results/Data/Model_dataset.csv")

    models_dict = get_models()


    first_class = ['Gain','Loss']
    second_class = ['Neutral']

    for c1 in first_class:
        for c2 in second_class:


            file_name = c1+' versus '+c2
            dataset_per_model, threshold = process_dataset_per_model(Dataset,c1,c2)

            dataset_per_model = dataset_per_model.drop(columns= ['Type','Arm'])

            dataset_per_model = null_to_median(dataset_per_model)

            X = dataset_per_model .loc[:, dataset_per_model .columns != 'Label']
            y = dataset_per_model ['Label']


            finding_hyperparameters(X,y,c1,c2,models_dict)


parameter_tunning()

