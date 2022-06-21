# Mixture Density Regression analysis of iHS on human populations

Pipeline used to test the association between the integrated haplotype Score (iHS; Voight *et al*., 2006) and multiple genomic factors by using a Mixture Density Regression (MDR) approach.

## Instructions

In order to run the MDR, you need to copy the script "MDR_script.py" and the folder "data". This folder include data for iHS and genomic factors for:

1. Five human populations from the [1000 genomes project](https://www.internationalgenome.org/):
	- Yoruba (YRID)
	- Toscani (TSID)
	- Utah residents (CEUD)
	- Han Chinese (CHBD)
	- Peruvians (CHBD)

2. Five window sizes:
	- 50 kb
	- 100 kb
	- 200 kb
	- 500 kb
	- 1000 kb

To run the script, you have to type the command `echo +x MDR_script.py` in the terminal and then you can use run it by typing `./MDR_script.py`. 

The script will automatically generate a "results" folder where figures and tables of each population and window size will be saved.

## References

1. Voight, B. F., Kudaravalli, S., Wen, X., & Pritchard, J. K. (2006). A map of recent positive selection in the human genome. PLoS biology, 4(3), e72. https://doi.org/10.1371/journal.pbio.0040072
