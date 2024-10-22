### Information  
The following repository accompanies the paper 
Lucic, B., Chen, H., Kuzman, M., Zorita. E. et al. Spatially clustered loci with multiple enhancers are frequent targets of HIV-1 integration. Nat Commun 10, 4059 (2019) [doi:10.1038/s41467-019-12046-3](https://www.nature.com/articles/s41467-019-12046-3)


They are supplementary material to the paper, and their usage is explained in the methods and supplementary methods of the paper.  

## Short description of the R scripts:  

Master.Rmd - the script that calls all the other scripts and contains explanations.  	 

1_RedefinitionOfRIGs.Rmd	 

2_AdditionalFile1.R	  

3_Supplementary figure 1a.R	  

4_Supplementary figure 1b.R	  

5_AllFigures.R	  

Brady_Integration_Sites.md	   

ROC_figures.R	  

Supplementary figure 3.R	  

## Short description of the files and data:  

IS.txt - the final list of integration sites used in the paper, with the datasets marked by authors who provided the data.	 

is.Robj	- integration sites used in the paper as an R object, this file is used in the scripts provided here and also contains integrations merged together in a different way, as is explained in the paper. The data set we provided in the paper (Lucic data set) was initially provided in two separate experiments, that can be found here under the names Lucic sorted and Lucic unsorted, but all of the integration sites are used together (as Lucic data set).  

IS_withaddedvalues_075.Robj	- integration sites and values of epigenetic features for the calculation of ROC analysis.   

RS_random_matched_controls.Robj	- random matched control sites and values of epigenetic features for the ROC analysis, I generated 10 random matched controls for every true site in the paper. See [Supplementary methods](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-12046-3/MediaObjects/41467_2019_12046_MOESM2_ESM.pdf) for details.     

genidf.Robj	- coordinates of genes used in the paper (hg19)   

roc.Robj - values for the ROC analysis plot    

Replicates.RDS - more random matched controls (100) generated for the review process. For details see [Peer review file](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-12046-3/MediaObjects/41467_2019_12046_MOESM2_ESM.pdf) .  	 
## List of R packages (dependencies for the .Robj and RDS files and scripts):   
R version >= 3.5  

## R libraries
```{r}
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install(c("GenomicRanges","Sushi","data.table","digest","stringr","pheatmap"), version = "3.8");
install.packages("gplots", repo="https://cran.cnr.berkeley.edu/")
install.packages("ggplot2", repo="https://cran.cnr.berkeley.edu/")
install.packages("RColorBrewer", repo="https://cran.cnr.berkeley.edu/")
install.packages("gridExtra", repo="https://cran.cnr.berkeley.edu/")
install.packages("scales", repo="https://cran.cnr.berkeley.edu/")
install.packages("rmarkdown", repo="https://cran.cnr.berkeley.edu/")
install.packages("tidyr", repo="https://cran.cnr.berkeley.edu/")
```

