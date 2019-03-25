# statistics
Tools for statistical analysis and testing.

## Introduction

This repository contains Python scripts I wrote that provide statistical tests and tools. The Jupyter notebook starting with a number for 1 to 13 correspond to an introductory course to statistics for biologists. Many areas are covered so the methods describes in these notebooks extend beyond the field of biology. These notebooks serve mostly the purpose of learning and teaching, so I mainly use my own algorithms to details the steps in statistical calculations. These algorithms work fine with datasets containing a few thousand points but may take too long to run on large datasets.

I will emphisize on a few notebooks:
* the notebook `4-linear_regression.ipynb` contains an rather complete study of linear regression, including verification of assumptions and a rather verbose output containing ample information to judge the quality of the linear regression.
* the notebook `6-anova.ipynb` presents the analysis of variance performed on data acquired using the completely randomised design. The following notebooks `7-randomised_block_design.ipynb` and `8-factorial_experiments.ipynb` extend ANOVA to other experimental designs. Again, results are wrapped in a verbose output to allow a extensive assessment of the statistical model.

## Description of files
"mannwhithney_exact_v2.py" performs an exact Mann-Whitney test (i.e. it does not use the normal distribution approximation) but it
does not handle ties so make sure your data does not have any).

"table_mw.txt" contains the U critical values and Wilcoxon cumulative values for different couples of samples sizes. Reading and
parsing is performed by "mannwhithney_exact_v2.py".

"generate_MW_table.py" was used to generate the data stored in "table_mw.txt". Computation time becomes very long for sample
sizes > 11.

"dunn_test.py" performs a Dunn's test (post-hoc non parametric test, can be used in conjonction with Kruskal-Wallis test to pinpoint which group is different from the others). The file was modified from a script found [here](https://gist.github.com/ricardoV94).
