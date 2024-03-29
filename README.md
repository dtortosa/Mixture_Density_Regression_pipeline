# Mixture Density Regression-based analysis of iHS in human populations

Pipeline used to test the association between the integrated haplotype Score (iHS; Voight *et al*., 2006) and multiple genomic factors across the human genome by using a Mixture Density Regression (MDR) model.

## Instructions

In order to run the MDR, you need to have the script [`MDR_script.py`](/MDR_script.py) and the folder [`data`](/data) in the same directory. The [`data`](/data) folder includes data about the average of iHS and the density (i.e., number) of iHS data points ([`mean_ihs_gene_windows`](/data/mean_ihs_gene_windows)), along with genomic factors ([`predictors_gene_windows`](/data/predictors_gene_windows)) across genomic windows centered around human coding genes. We consider five window sizes: 50 kb, 100 kb, 200 kb, 500 kb and 1,000 kb. 

iHS was calculated across five human populations from the [1000 genomes project](https://www.internationalgenome.org/):

- Yoruba (YRI)
- Toscani (TSI)
- Utah residents (CEU)
- Han Chinese (CHB)
- Peruvians (PEL)

In order to estimate the correlation between iHS and the genomic factors, you have to type the command `chmod +x MDR_script.py` in the terminal and then you can run the script by typing `./MDR_script.py > MDR_script.out`. Note that you need python 3.8 and the corresponding python libraries installed in order to run this script (you can find them in the Imports section of [`MDR_script.py`](/MDR_script.py)). 

The script will automatically generate a [`results`](/results) folder where figures and tables of each population and window size will be saved (5 populations x 5 windows sizes = 25 tables and 25 figures). You can check that all files have been created in [`MDR_script.out`](/MDR_script.out) once the script has been run. The figures will be saved as png, while the tables will be saved as csv.

See [Salazar-Tortosa *et al*. (2023)](https://doi.org/10.1093/gbe/evad170) for further details.


## References

1. Voight, B. F., Kudaravalli, S., Wen, X., & Pritchard, J. K. (2006). A map of recent positive selection in the human genome. PLoS biology, 4(3), e72. https://doi.org/10.1371/journal.pbio.0040072
2. Salazar-Tortosa, D. F., Huang, Y-F., Enard, D. (2023). Assessing the presence of recent adaptation in the human genome with Mixture Density Regression. Genome Biology and Evolution, evad170. https://doi.org/10.1093/gbe/evad170

