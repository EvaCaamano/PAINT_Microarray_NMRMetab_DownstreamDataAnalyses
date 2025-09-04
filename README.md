# PAINT_Microarray_NMRMetab_DownstreamDataAnalyses
Repository with main statistical analyses for the project PAINT - Analyses led by the Computational Biology Facility at the Univerisity of Liverpool. Samples and Study led by the Liverpool Women's Hospital Neonatal Team.

Omics datasets were produced by the Centre for Genomic Research at the Univeristy of Liverpool (microarray data) and C. Slupsky's team at UCL (NMR Metabolomics data). This repository only includes analyses undertaken by the Computational Biology Facility, downstream from processing.The study followed neonatal patients in two data points day 3 and day 10 of life.

File AnalyticalSteps contain the main steps undertaken for the analysis in each omics dataset:
1. Variance partition to understand the variance explained by the two timepoints, patient ID, Gestation time at birth and Sex of the patients
2. Differential analysis - via limma using a dupplicated correlations approach to account for patient variability. Three models were fitted exploring Time only, adding Sex as covariate and adding gestation at birth and sex as covariates.
3. Further visualisations - e.g. PCA, Heatmap if appropriate
4. For microarray results functional enrichment analysis with package fgsea on Reactome database
