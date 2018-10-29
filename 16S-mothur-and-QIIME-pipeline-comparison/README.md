# mothur_16S
Pipeline for 16S analysis with mothur

This pipeline is based on mothur 16S MiSeq SOP. It is prepared for performing sample paralelization, and thus processes samples per separate. Script:
- **MOTHUR_run**: Includes the mothur pipeline. It also contains an optional trimmommatic pre-processing step. Inside this script mothur parameters can be customized and adapted depending on user's dataset. THis script can be executed for one sample or using parallelization.
