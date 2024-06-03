# nfds-tutorial
Simple tutorial of technical workflow to study negative frequency dependent selection (NFDS) in the accessory genome of Streptococcus pneumoniae. Covers the bioinformatic steps to go from shotgun whole genome sequencing raw reads to accessory gene frequencies. The accessory gene frequencies can be used to predict the post-vaccination population structure of Streptococcus pneumoniae. The implementation of this NFDS-based quadratic programming predictive model is also discussed.

## Table of Contents
1. [Introduction](#introduction)
2. [Data pre-processing](#data-pre-processing)
3. [Pangenome Analysis](#pangenome-analysis)
4. [Strain Typing](#strain-typing)
5. [Quadratic Programming Model](#quadratic-programming-model)
6. [Conclusion](#conclusion)

## Introduction

This tutorial relies heavily on conda environments so a conda installation is assumed. For quick and simple instructions of how to locally install conda with miniconda follow [this link](https://docs.anaconda.com/free/miniconda/#quick-command-line-install).

I will focus on For the steps in which the same command is generally run for multiple files (cleaning, annotations, assemblies, etc.) I use bash scripts to perform the step in batch, which is more useful and closer to how we run these analyses. They are especially helpful if you want to submit jobs on a computing cluster. These scripts are found in the [scripts](/scripts) folder within this repository.

## Data pre-processing

In this step we will go from raw reads to assemblies, doing appropiate quality control for each step. 

### Setting up main environment

Due to incompatibilities, two environments are used for this step: qc_assembly.yml and quast.yml. The qc_assembly environment is used for all cleaning and assembly steps, while quast is only used for the last step of assembly qc. Both of these files are located in the [envs](/envs) folder within this repository.

To create the conda environments for the read cleaning and assembly steps:

```bash
conda env create --file envs/qc_assembly.yml
```

then activate the environment:

```bash
conda activate qc_assembly
```

Now you should be able to run the data pre-processing steps. For these examples we will be using sample ERR065307 as a case study.

### Read cleaning with [Trimmomatic](https://github.com/timflutre/trimmomatic)

First we navigate to the folder that contains the example raw reads:

```bash
cd  ~/exercises/data_preprocessing/reads
```
Now, we can run the trimmomatic command:

```bash
trimmomatic PE -phred33 \
ERR065307_1.fastq.gz ERR065307_2.fastq.gz \
~/exercises/data_preprocessing/trimmed_reads/ERR065307_R1_trimmed.fastq.gz /dev/null \
~exercises/data_preprocessing/trimmed_reads/ERR065307_R2_trimmed.fastq.gz /dev/null \
SLIDINGWINDOW:4:20 MINLEN:30 ILLUMINACLIP:NexteraPE-PE.fa:2:30:10
```
The clean reads will be found in the 'trimmed_reads' folder and we can proceed to QC the new read files.

### Read quality control with [FastQC](https://github.com/s-andrews/FastQC)

Similarly, we navigate to the directory that contains the clean reads:

```bash
cd  ~/exercises/data_preprocessing/trimmed_reads
```
Now, we can run the FastQC command:

```bash
fastqc -o ~/exercises/data_preprocessing/fast_qc ERR065307_R1_trimmed.fastq.gz ERR065307_R2_trimmed.fastq.gz
```
And to summarize all FastQC reports into one, we run MultiQC:

```bash
multiqc ~/exercises/data_preprocessing/fast_qc
```
Now you can assess the read quality accross all samples.

### Assembly with [Unicycler](https://github.com/rrwick/Unicycler?tab=readme-ov-file#quick-usage)

Assemblies are one of the most computationally intensive steps in the pipeline. Thus, they are usually run in high performance computing clusters. From the folder containing the clean reads you can run Unicycler:

```bash
unicycler -1 ERR065307_R1_trimmed.fastq.gz -2 ERR065307_R2_trimmed.fastq.gz -o ~/exercises/data_preprocessing/assemblies/ERR065307_assembly
```
Unicycler generates many output files, and there will be one folder per sample. The resulting assembly will be named assembly.fasta

### Assembly quality control with [Quast](https://github.com/ablab/quast)

## Pangenome analysis

### Annotate assemblies with [Prokka](https://github.com/tseemann/prokka)

### Get pangenome gene definitions with [Roary](https://github.com/sanger-pathogens/Roary)

### Clean Roary gene definitions with CLARC (obtain final presence absence matrix for accessory genes)

## Strain typing 

### Strain typing using [PopPUNK](https://github.com/bacpop/PopPUNK)

## Quadratic programming model 

### Predicting post-vaccine population structure with quadratic programming (QP) model

## Conclusion


