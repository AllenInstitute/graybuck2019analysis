## graybuck2019analysis
Scripts used for data analysis for Graybuck, et al. (2019).

These scripts are currently available as-is, including references to local files and directories that may not be useful to external users.

## Alignment
An example of alignment scripts for a single 50bp paired-end MiSeq run is stored in the `alignment_example/` directory.  

Parameters for the locations of alignment tools and genomes must be provided in `00_config.sh`.

The full alignment procedure as written requires the following tools:  
* bowtie: http://bowtie-bio.sourceforge.net/index.shtml  
* cutadapt: https://github.com/marcelm/cutadapt/  
* trimgalore: https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/  
* samtools: http://www.htslib.org/  
* preseq: http://smithlabresearch.org/software/preseq/  
* bamToBed: Part of the UCSC kent tools collection http://hgdownload.soe.ucsc.edu/admin/exe/  
* Picard Tools: https://broadinstitute.github.io/picard/  

## scATAC-seq analysis

The primary scATAC-seq analysis pipeline is written in R, and is stored in the `gb_analysis` directory.  

This analysis pipeline requires many R packages, which can be checked and installed with:
```
gb_cran_reqs <- c("dplyr",
                  "purrr",
                  "ggplot2",
                  "ggbeeswarm",
                  "dbscan",
                  "Rtsne",
                  "data.table",
                  "feather",
                  "Matrix",
                  "matrixStats",
                  "WGCNA",
                  "rbamtools",
                  "BiocManager",
                  "devtools")

missing_cran_reqs <- setdiff(gb_cran_reqs, installed.packages())
install.packages(missing_cran_reqs)

gb_bioc_reqs <- c("GenomicAlignments",
                  "GenomicRanges",
                  "rtracklayer",
                  "limma",
                  "Rsamtools")

missing_bioc_reqs <- setdiff(bg_bioc_reqs, installed.packages())
BiocManager::install(missing_bioc_reqs)

gb_ghub_reqs <- c("AllenInstitute/scrattch.hicat",
                  "AllenInstitute/scrattch.vis",
                  "AllenInstitute/lowcat",
                  "JinmiaoChenLab/Rphenograph")

missing_ghub_reqs <- setdiff(sub(".+/","",gb_ghub_reqs), installed.packages())
devtools::install_github(missing_ghub_reqs)
```
Separate scripts for joint analysis of data from Graybuck, et al. (2019) and Prefrontal Cortex data from Cusanovich and Hill, et al. (2018) are stored in `supp_gb_cu_joint_analysis/`.

## Level of Support

We are not currently supporting this code, but simply releasing it to the community AS IS but are not able to provide any guarantees of support. The community is welcome to submit issues, but you should not expect an active response.

## License

This code is released under the Allen Institute Software License. See file LICENSE for details.

## Contributions

Any contributions to this repository are subject to the Allen Institute Contribution Agreement. See file CONTRIBUTING.md for details.
