from compare_models import *
from shap_analysis import *
from dataset_process_main import *
from CCL_main import *
from correlation_analysis import *
from paralogs_analysis import *
def main():

    '''
    processing and constructing dataset for ML model
    '''
    dataset_processing_main()

    '''
    Assessing the performance of 5 different ML models, for all the models including:
    Gain model:  Gain versus Neutral
    Loss model: Loss versus Neutral
    Gain model: Gain versus Rest (i.e. Loss and Neutral)
    Loss model: Loss versus Rest (I.e. Gain and Neutral)
    Results are saved as Figure S3
    '''
    compare_ML_methods_performance()



    '''
    Using SHAP algorithm to reveal contributing features per each model
    '''
    SHAP_analysis_of_features_contributing()


    '''
    Construct CCL dataset for validation
    '''
    CCL_calculate()

    '''
    Calculate correlation between features and arm aneuploidy
    '''
    calculate_correlations()


    '''
    Paralog analysis
    '''
    analyze_lost_paralogs()
main()