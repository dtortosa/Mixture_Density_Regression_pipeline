# Mixture_Density_Regression_pipeline

Pipeline used to test the association between the selection statistic iHS and multiple genomic factors by using a Mixture Density Regression (MDR) approach.

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

To run the script, you have to type the commmand `echo +x MDR_script.py` and then you can use run it as `./MDR_script.py`. 

The script will automatically generate a "results" folder where figures and tables of each population and window size.
