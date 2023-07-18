# ML_Aneuploidy_interpretation
This repository contains data, code, and analysis as described in the paper "Machine-learning analysis of factors that shape cancer aneuploidy landscapes reveals an important role for negative selection" by Juman Jubran, Rachel Slutsky, Nir Rozenblum, Lior Rokach, Uri Ben-David and Esti Yeger-Lotem

## Table of contents
* [Article results](#article-results) 
    * [Pipeline code](#pipeline-code)
    * [Plots](#plots)
* Data availability

### Article results
Results are found in their respective "article_results/Figure X" folder.

#### Pipeline code
the main code was written in python. It includes data processing, machine learning construction and validation. 
The pipeline code can be reviewed and used by running main.py from src folder.
(The dataset and the code is being updated and uploaded)

#### Plots
Plots were generated via python as part of the pipeline and via R.

| Figure   | Panel | Script                                         |
|----------|-------|------------------------------------------------|
| Figure1  | A     | build_aneuploidy_table_GISTIC.py (pipeline)    |
| Figure2  | A-D   | shap_analysis.py (pipline)                     |
| Figure4  | B     | Fig4B_13q_gain_prevalence_barplot              |
|          | D,F   | violin_plots.R                                 |
| Figure5  | D,E   | violin_plots.R                                 |
|          | G,H   | Fig5G-H_cooccurrence.R                         |
| FigureS3 | A-D   | curves - compare_models.py (pipline)           |
|          |       | boxplot - FigS3_models_performance.R           |
| FigureS4 | A,B   | shap_analysis.py (pipeline)                    |
| FigureS6 | A-D   | shap_analysis.py (pipeline)                    |
