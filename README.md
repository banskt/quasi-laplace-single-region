# Simulation pipeline for comparing B-LORE against full MCMC in a combined locus in a combined study.

B-LORE uses quasi-Laplace approximation for multiple logistic regression.
It benefits over existing MCMC logistic regression from two aspects: 
 1. combining information over multiple loci
 2. combining information over multiple studies

In the limit of single locus and single study, MCMC logistic regression should ideally perform as good as B-LORE.
We explore this limit in this simulation pipeline. We run B-LORE on 25 loci over 5 studies.
We combine them to a single big locus in a single big study,
and compare the quasi-Laplace multiple logistic regression with other methods:
 * multiple probit regression (piMASS, GEMMA) using point-normal prior
 * multiple logistic regression with linear approximation (FINEMAP) using point-normal prior
 * multiple linear regression (piMASS, GEMMA) using point-normal prior
 * simple logistic regression (SNPTEST)

In our simulations, we vary the following parameters:
 * heritability (h)
 * maximum number of causal SNPs allowed (c)
 * ratio of cases to controls (&lambda;)

We compare the posterior inclusion probabilities (PIPs) with the following settings:
 - [ ] varying h = 0.2, 0.4, 0.6 and 0.8 (fixed c = 2, &lambda; = 1)
 - [ ] varying c = 2, 3, 4 and 5 (fixed h = 0.4, &lambda; = 1)
 - [ ] varying &lambda; = 0.25, 0.50, 0.75, 1.00 (fixed h = 0.4, c = 2)
 - [ ] including covariates (fixed h = 0.4, c = 2, &lambda; = 1)

## Method
We used the original genotype and sample files from the five GerMIFS studies (total 13082 patients)
to select 5000 SNPs distributed over 25 loci as input.
We performed the following tasks in this pipeline:
 * create the loci using SNPs which are common to all studies
 * combine the loci into a single big locus for all studies
 * combine the studies to a single big study
 * perform SNPTEST / META analysis with the original phenotype
 * find the LD matrix for each loci (requires META analysis for proper ordering of SNPs)
 * simulate the phenotype using different parameters
 * run SNPTEST / META, FINEMAP, PIMASS, GEMMA and BLORE on the data
 * plot the precision and recall values for each setting

## Input files
The pipeline expects input genotype and phenotype to be organized in the following way:
```
{DOSAGEDIR}
   -- {STUDY}
         -- {LOCUSPREFIX}.gen
         -- {LOCUSPREFIX}.map
         -- {STUDY}{SAMPLEPREFIX}.sample
```
There could be multiple ```{STUDY}``` folders in the ```{DOSAGEDIR}``` and multiple ```{LOCUSPREFIX}``` within each study. 
However the ```{LOCUSPREFIX}``` should be exactly same in each study. 

The genotype and sample files are in the Oxford format.
Additionally, the sample file should contain covariates called 'sex' and 'age'. 
The binary phenotype should be in the last column, which should be named 'pheno'.

## How to run
The simulation can be run from the ```pipeline``` folder. Update the ```CONFIG``` and ```PATHS``` file and run:
```
./01_create_loci.sh
./02_find_loci_with_max_snps.sh
./03_make_single_biglocus.sh
./04_run_snptest_meta.sh
./05_create_ldmatrix.sh
./06_simulate.sh CONFIG
./07_getplots.sh CONFIG
```
