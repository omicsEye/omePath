# deepath
deepath is a generic tool for omics pathway enrichment analysis

* **Citation:** if you use the deepath software, please cite our manuscript: Rahnavard et al. **deepath: a generic tool for omics enrichment analysis**.
* This page provides a quick toturial (workshop orinted) infomration to start and use *deepath*. 


If you have questions, please submit it as an issue using [deepath issue tracker](https://github.com/omicsEye/deepath/issues).
--------------------------------------------

## Contents ##
* [Description](#description)
* [Requirements](#requirements)
* [Installation](#installation)
* [How to Run](#how-to-run)
    * [Input Files](#input-files)
    * [Output Files](#output-files)
    * [Run a Demo](#run-a-demo)
    * [Options](#options)
* [Visualization](#visualization)
* [Troubleshooting](#troubleshooting)

## Description ##

deepath was developed to ... 

## Requirements ##

deepath is an R package that can be run on the command line or as an R function. It requires the following R packages included in Bioconductor and CRAN (Comprehensive R Archive Network). Please install these packages before running deepath.

* Only for Windows OS:
    * ``Install Rtools:`` https://cran.r-project.org/bin/windows/Rtools/
    * We tested using `Rtools35.exe` version, and in the future we will ensure `deepath` works with the latest recommended version.

* CRAN packages
    * [dplyr: A Grammer of Data Manipulation](https://cran.rstudio.com/web/packages/dplyr/index.html)
    * [ggplot2: Create Elegant Data Visualizations Using the Grammer of Graphics] https://cran.rstudio.com/web/packages/ggplot2/index.html)
    * [readr: provide a fast and friendly way to read rectangular data (like 'csv', 'tsv', and 'fwf')](https://cran.r-project.org/web/packages/readr/index.html)
    * [downloader: Download Files over HTTP and HTTPS](https://cran.r-project.org/web/packages/downloader/index.html)
    * [logging: R logging package](https://cran.rstudio.com/web/packages/logging/index.html)
    * [data.table: Fast aggregation of large data](https://cran.rstudio.com/web/packages/data.table/index.html)
    * These packages can be installed in R with ``install.packages('dplyr')`` or from the command line ``$ R -q -e "install.packages('dplyr', repos='http://cran.r-project.org')"`` individually (for those packages which you do not yet have installed) or as a set by providing the complete list as a vector.

## Installation ##

deepath can be installed as an R package and run as an R function. You will need to install the deepath dependencies.

### From R ###


### Install deepath in RStudio ###

1. Install devtools : 
    * ``> install.packages('devtools')``
    * ``>library(devtools)``
2. Install the Bioconuctor dependencies: 
    * ``> install.packages('BiocManager'); library('BiocManager');``
    * ``> BiocManager::install('limma')``
3. Install the CRAN dependencies:
    * ``> install.packages(c('future', 'downloader', 'reader', 'backports', 'gsEasy','pscl','pbapply','car','nlme','dplyr','vegan','chemometrics','ggplot2','pheatmap','cplm','hash','logging','data.table'), repos='http://cran.r-project.org')``
4. Install deepath (and also all dependencies from CRAN): 
    * ``> devtools::install_github('omicsEye/deepath', force = TRUE)``

### Download the mapping database ###

Users can bring their own mapping files (pathways-omics) with the following format:

* One of the columns should be labled as `Pathway`. ....
* One of the columns should be labled as `Feature`.
* The table should be a  tab-delimited text file (`tsv`).
* Input file columns should have consistent labels.
    
We provide mapping files (pathways-feature) for four main omics:

* Metagenomics: Microbial genes pathways (GO terms - UniRef90)

* The Molecular Signatures Database (MSigDB), C5: ontology gene sets including 14765 gene sets
    * Please download the preprocessed mapping file at [HERE](https://www.dropbox.com/s/zlobyzn92r43nqy/c5.all.v7.1.entrez.tsv?dl=0)


* Metabolomics: (Metabolic Pathways - Metabolites HMDBID)
    * Please download the preprocessed mapping file at [HERE](https://www.dropbox.com/s/rzx6hsq7xmr87te/smpdb_metabolites.tsv?dl=0)
    

* Proteomics: coming soon


### Input files format ###

* all input files should be tab-delimited formatted  
* deepath requires a mapper file that can be downloaded from [previouse step](#download_the_mapping_database) an input files.

#### Two approaches of run *deepath* ####


1. Providing the score file and mapping file:

In this approach, a tab-delimited text file or a R data frame with row names being the features need to be provided.
The file should have a column which will be used as the score for enrichment analysis. 

[deepath demo](https://github.com/omicsEye/deepath/tree/master/demo)

|               |  logFC         |  statistic  |  P.Value
----------------|----------------|-------------|-------------
HMDB00696       |  -0.066838102  |  60         |  0.513722581
HMDB00191       |  -0.167002899  |  68         |  0.84283599
HMDB00148       |  -0.106438214  |  60         |  0.513722581
HMDB00168       |  -0.029952907  |  69         |  0.887385935
HMDB00641       |  -0.040854959  |  64         |  0.670659533
HMDB00177       |  -0.050286514  |  55         |  0.347357919
HMDB00517       |  -0.0184334    |  70         |  0.932300503
HMDB00182       |  -0.001386226  |  70         |  0.932300503
HMDB00883       |  0.04519753    |  74         |  0.932300503
HMDB00687       |  0.049785325   |  73         |  0.977402191
HMDB00172       |  0.04173989    |  74         |  0.932300503
HMDB00159       |  -0.034784137  |  72         |  1
HMDB00158       |  0.034486909   |  76         |  0.84283599
HMDB00929       |  0.071698012   |  80         |  0.670659533
HMDB00162       |  -0.030734138  |  67         |  0.798745339
HMDB00725       |  0.070803166   |  92         |  0.265669584

and continues ... 

2. Providing a data and a metadata file to calculate a score for feature (currenlty logFC). providing data, metadata, and mapping file: this allows to score features using provides function by deepath such as logFC.

* Data (or features) file : This file is tab-delimited formatted with features as columns and samples as rows (the transpose is also okay).

* Metadata file: This file is tab-delimited formatted with metadata as columns and samples as rows (the transpose is also okay).
    
* Mapping files: maps pathways to memebrs (features).

NOTE: If running deepath as a function, the data inputs can be of type ``data.frame`` instead of a path to a file.


## Run ##

```
# load the library
library(deepath)


# call the function
deepath_results <- deepath(input_data,
                    output,
                    input_metadata = NA,
                    meta = NA,
                    case_label = NA,
                    control_label = NA,
                    score_col = 'logFC',
                    pval_threshold = 0.05,
                    fdr_threshold = NA,
                    Pathway.Subject = NA,
                    method = 'ks',
                    min_member = 2,
                    mapper_file = NA,
                    do_plot = TRUE,
                    pathway_col = "Pathway",
                    feature_col = "Feature")
                    

```

Required option for both running approaches:

* `mapper_file` should have at least two columns to map pathways ("Pathway" column by default) to feature ("Feature" column by default). Different mapper files could use different column labels for pathways (e.g, 'Pathway', 'KEGG', and 'GO') and 
feature (e.g, 'Feature', 'metabolites', and 'HMBD_ID', 'Gene', 'Ensembl_ID')

* `pathway_col` identifies pathways column in the `mapper_file`, default is "Pathway"

* `feature_col` identifies features column in the `mapper_file`, default is "Feature"


In running approach 1 (providing score file and mapping file), following data files need to be specified:

* `input_data` will be your score file with row names your features and a column for your scores (e.g. 'logFC' and 'coef')

* `score_col` identifies your score column in th e`input_file`, default is `logFC`


In running approach 2 following data files need to be specified::

* `input_data` will be omics profile data with N rows by D columns where N are samples and D are features. 

* `input_metadata` identifies your metadata file with N rows by M columns where N are samples and M are metadata. 


* `meta` is the name of teh column in metadata that has the case and control groups and need to 
be used for calculating score (e.g., Diagnosis).

* `case_label` is the label used for the case group in `meta` column (e.g., 'UC' or 'Ulcerative Colitis') 

* `control_label`is the label used for the control group in `meta` column (e.g., 'non-IBD') 


### Output Files ###

deepath generates two types of output files: data and visualization.

1. Data output file
    * ``enrichment_stats.tsv`` : This file contains all of the enrichment analysis results ordered by increasing q-value.

pathway                          |  pval         |  fdr          |  n   |  set_enrichment_score
---------------------------------|---------------|---------------|------|----------------------
Arginine and Proline Metabolism  |  0.333333333  |  0.333333333  |  10  |  0.425686617
Histidine Metabolism             |  0.179104478  |  0.268656716  |  10  |  0.502342185
Purine Metabolism                |  0.001369986  |  0.004109959  |  10  |  0.701651617


    * ``deepath.log`` : This file contains all of the debug information for the run. It includes all settings, warnings, errors, and steps run.
2. Visualization output files
    * ``enrichment_plots.pdf`` : This file contains a heatmap of the significant associations.

### Run a Demo ###

Example input files can be found in the tests folder of the deepath source. 



### Options ###

Run deepath help to print a list of the options and the default settings.




## Visualization ##

`deepath` generates two enrichment plot per pathways (rank based and score based) in deepath that visualize the outputs and provide ggplot2 plots that can be used to generate manuscript/report quality figures. the m2 path returns a result variable which contains three variables 1) stats table 2) a list or ranked based enrichment plots, and 3) a list of score based enrichment plots.
* ``deepath ``: this function generates an enrichment plot of all pathways reported by deepath and have the following parameters:

``output_path``: the path to the deepath output.

``do_plot``: a parameter to default is TRUE to generate the plots

``pval_threshold``: a threshold used to visualize only pathways with enrichment p-value less or equal to it, and the default is 0.05.

``fdr_threshold``: a threshold used to visualize only pathways with enrichment q-value (FDR) less or equal to it, and the default is NA. If the value is not NA then overwrites the ``pval_threshold`` condition.

[[/img/enrichment_rank_based_plot.png|alt=Enrichment Plot]]

## Troubleshooting ##

1. Question: When I try to install the R package I see errors about dependencies not being installed. Why is this?
    * Answer: Installing the R package will not automatically install the packages deepath requires. Please install the dependencies and then install the deepath R package.

2. Question: When I run as a function I see the error ``Error in library(deepath): there is no package called 'deepath'``. How do I fix this? 
    * Answer: Install the R package and then try loading the library again.


