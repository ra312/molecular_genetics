# molecular_genetics
We seek for patterns in genomic data from clinical trials.
A python script parse.py transforms a XML-file igv_session.xml into a Python object.
The main programme loci_gene.c reads the main pool of genes (all_test.txt) and the malignant genes (selected.txt).
Then we conduct statistical analysis (computing mode value, mean value, min and max values, dispersion of the gene values in the selected regions) of the drug induced values in the genes and compare the genomic regions according to the mean value and the dispersion.

