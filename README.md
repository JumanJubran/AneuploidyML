# ML_Aneuploidy_interpretation
This repository contains data, code, and analysis as described in the paper "Machine-learning analysis of factors that shape cancer aneuploidy landscapes reveals an important role for negative selection" by Juman Jubran, Rachel Slutsky, Nir Rozenblum, Lior Rokach, Uri Ben-David and Esti Yeger-Lotem

### The repository is under construction

## Table of contents
* [Article results](#article-results) 
    * [Pipeline code](#pipeline-code)
    * [Plots](#plots)
* [Data availability](#Data-availability)

### Article results
Results are found in their respective "article_results/Figure X" folder.

#### Pipeline code
the main code was written in python. It includes data processing, machine learning construction and validation. 
The pipeline code can be reviewed and used by running main.py from src folder.

#### Plots
Data to generate each panel is presented in "Data source file.xlsx" and in their respective "article_results/Figure X" folder.
Plots were generated via python as part of the pipeline and via R.

| Figure    | Panel | Script                                                                                                                   |
|---------- |-------|------------------------------------------------------------------------------------------------------------------------- |
| Figure1   | A     | build_aneuploidy_table_GISTIC.py (pipeline)                                                                              |
| Figure2   | A-D   | shap_analysis.py (pipline)                                                                                               |
| Figure2   | E     | correlation_analysis.py (pipline)                                                                                        |
| Figure3   | A-E   | CCL_main.py (pipline)                                                                                                    |
| Figure3   | F     | correlation_analysis.py (pipline)                                                                                        |
| Figure4   | B     | Fig4B_13q_gain_prevalence_barplot                                                                                        |
|           | D,F   | violin_plots.R                                                                                                           |
|           | H,J   | Fig4H-J_boxplots                                                                                                         |
| Figure5   | A     | correlation_analysis.py (pipline)                                                                                        |
| Figure5   | D,E   | violin_plots.R                                                                                                           |
|           | G,H   | Fig5G-H_cooccurrence.R                                                                                                   |
| FigureS3  | A-C   | Correlation analysis: 'article_results\Figure S3\correlation analysis.py                                                 |
|           |       | Plots: 'article_results\Figure S3\Fig_S3_plots.r                                                                         |
| FigureS4  | A-B   | curves - compare_models.py (pipline)                                                                                     |
|           |       | boxplot - FigS3_models_performance.R                                                                                     |
| FigureS5  | A,B   | shap_analysis.py (pipeline)                                                                                              |
| FigureS6  |       | correlation_analysis.py (pipline)                                                                                        |   
| FigureS7  |       | construct dataset with other threshold: src\other_threshold_percentage.py                                                |  
|           |       | To find best model and to plot SHAP, run compare_models.py and SHAP analysis respectively on the relevant dataset.       | 
| FigureS8  | A     | parameter_tunning.py (pipline)                                                                                           |
| FigureS9  | A     | curves - compare_models.py (pipline)                                                                                     |
| FigureS9  | B-C   | shap_analysis.py (pipeline)                                                                                              |
| FigureS17 | A-B   | paralogs_analysis.py (pipline)                                                                                           |

#### Data availability
Dataset folder is required to build aneuploidy dataset to construct the ML model.
The dataset is available on Zenode - https://zenodo.org/record/8199048
