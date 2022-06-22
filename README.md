# Mixture Density Regression analysis of iHS on human populations

Pipeline used to test the association between the integrated haplotype Score (iHS; Voight *et al*., 2006) and multiple genomic factors across the human genome by using a Mixture Density Regression (MDR) approach.

## Instructions

In order to run the MDR, you need to have the script `MDR_script.py` and the folder `data` in the same directory. The `data` folder includes data of average iHS, the density (i.e., number) of iHS data points and genomic factors across windows centered around human coding genes. We consider five window size: 50 kb, 100 kb, 200 kb, 500 kb and 1,000 kb. 

iHS was calculated across five human populations from the [1000 genomes project](https://www.internationalgenome.org/):

- Yoruba (YRID)
- Toscani (TSID)
- Utah residents (CEUD)
- Han Chinese (CHBD)
- Peruvians (PELD)

To run the script, you have to type the command `echo +x MDR_script.py` in the terminal and then you can use run it by typing `./MDR_script.py`. Note that you need python 3.8 and the corresponding python libraries installed in order to run this script. 

The script will automatically generate a `results` folder where figures and tables of each population and window size will be saved. The figures will be saved as pdf, while the tables will be compressed as text files, being easy to open in a spreadsheet.


## References

1. Voight, B. F., Kudaravalli, S., Wen, X., & Pritchard, J. K. (2006). A map of recent positive selection in the human genome. PLoS biology, 4(3), e72. https://doi.org/10.1371/journal.pbio.0040072
