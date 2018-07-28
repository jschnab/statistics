# statistics
Tools for statistical analysis and testing.

This repository contains Python scripts I wrote that provide statistical tests and tools.

"mannwhithney_exact.py" performs an exact Mann-Whitney test (i.e. it does not use the normal distribution approximation) but it
does not handle ties so make sure your data does not have any.

"table_mw.txt" contains the U critical values and Wilcoxon cumulative values for different couples of samples sizes. Reading and
parsing is performed by "mannwhithney_exact.py".

"generate_MW_table.py" was used to generate the data stored in "table_mw.txt". Computation time becomes very long for sample
sizes > 11.
