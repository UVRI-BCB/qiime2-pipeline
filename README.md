# ![galaxyuvri-ea/qiime2-pipeline](/assets/uvrilogo.png)

# 16S rRNA gene microbiome data analysis pipeline using qiime2, implemented in Nextflow

A qiime2-based workflow using the Nextflow workflow manager.

## Dependencies

* fastQC v0.11.8 (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 
* multiQC v1.0.dev (https://multiqc.info/)
* qiime2 v2018.11 (https://docs.qiime2.org/2018.11/) 
* picrust2 v2.1.4 (https://github.com/picrust/picrust2)
* Nextflow v18.10.1 (https://www.nextflow.io/)

## Installation

Clone the pipeline repository in a desired directory:

```
git clone https://github.com/galaxyuvri-ea/qiime2-pipeline.git
```

The repository includes:

* the Nextflow script, main.nf,
* the configuration files for nextflow, nextflow.config
* a folder (conf) containing config files to specify computing infranstructure parameters

**Note**: the config files include the parameters that have been used in this exercise, one needs to appropriately change these most especially parameters pointing to file locations.

## Basic usage:

```
nextflow run qiime2-pipeline -profile base -reads 'path/to/reads' -reference 'path/to/reference' -man 'path/to/manifest file' -denoiser 'dada2' 
```
Where ;
* `reads` are path to paired-end reads, 
* `reference` path to pre-trained qiime2 classifier
* `man` path to manifest file to be used for importing data into qiime2 artifact
* `denoiser` string specifying denoiser to be for infering amplicon sequence variants

The pipeline can also be run without installing it first. In addition to supplying the mandatory arguments, we need to specify the repository name and everything else is taken care of by nextflow.

```
nextflow run galaxyuvri-ea/qiime2-pipeline -profile base -reads 'path/to/reads' -reference 'path/to/reference' -man 'path/to/manifest file' -denoiser 'dada2' 
```

## Singularity

In order to run the pipeline with singularity, we need to specify this using `-profile singularity`. Note that the singularity config may be changed to include specifications for resources assigned to individual processes. This needs to be appropriate to the computing facility used for analysis.

```
nextflow run qiime2-pipeline -profile singularity -reads 'path/to/reads' -reference 'path/to/reference' -man 'path/to/manifest file' -denoiser 'dada2' 
```
