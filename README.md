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

## Data pre-processing

In this step we will go from raw reads to assemblies, doing appropiate quality control for each step. 

### Setting up environment

### Read cleaning with [Trimmomatic](https://github.com/timflutre/trimmomatic)

### Read quality control with [FastQC](https://github.com/s-andrews/FastQC)

### Assembly with [Unicycler](https://github.com/rrwick/Unicycler?tab=readme-ov-file#quick-usage)

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


