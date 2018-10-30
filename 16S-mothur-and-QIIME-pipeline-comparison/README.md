# Comparison of mothur and QIIME for the analysis of rumen microbiota composition based on 16S rRNA amplicon sequences.
Pipelines for 16S analysis with mothur (v. 1.39.5) and QIIME (v. 1.9.1). Reference database can be chosen between GreenGenes (13_8) and SILVA (rel. 132). 

- **MOTHUR_run**: This script contains the mothur pipeline based on mothur [16S MiSeq SOP](https://www.mothur.org/wiki/MiSeq_SOP). It is prepared for performing sample paralelization, and thus processes samples per separate. It also contains an optional trimmommatic pre-processing step. Inside this script mothur parameters can be customized and adapted depending on user's dataset.
- **QIIME_run**: This script contains the QIIME pipeline for closed reference otu picking approach. It also contains an optional trimmommatic pre-processing step. As for mothur pipeline, QIIME parameters can be customized depending on user's dataset.
