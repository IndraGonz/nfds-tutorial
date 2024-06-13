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

I will focus on the most simple execution of each command. However, for the steps in which the same command is generally run for multiple files (cleaning, annotations, assemblies, etc.) I use bash scripts to perform the step in batch, which is more useful and closer to how we run these analyses. They are especially helpful if you want to submit jobs on a computing cluster. These scripts are found in the [scripts](/scripts) folder within this repository.

If you wish to run these exercises locally (or from a cluster), you can copy the repository:

```bash
git clone https://github.com/IndraGonz/nfds-tutorial.git
```
This tutorial is structured to be run from the repository, but of course I'm not the boss of you and you can do as you wish!

Finally, GitHub cannot host large files, so the raw input files used in this tutorial can be found in google drive. I will specify what files should be downloaded at each step and in what folder of the repository they should go in, to smoothly continue the tutorial.

With all that said, let's get started!

## Data pre-processing

In this step we will go from raw reads to assemblies, doing appropiate quality control for each step. 

### Setting up main environment

Due to incompatibilities, three environments are used for this step: `qc_assembly.yml, unicycler.yml, quast.yml`. The qc_assembly environment is used for all cleaning steps, unicycler is used for assembly and quast is only used for the last step of assembly qc. All of these files are located in the [envs](/envs) folder within this repository.

**Note:** We will be creating a lot of different environments for this tutorial. Like A LOT of environments so, bear with me. It's one way I've found to minimize things conflicting and crashing.

Anyways, to create the environments we use the conda create command, which must be run from the envs folder that contains all the precious .yml files. So, the unspoken pre-command for all the 'conda create...' commands in this tutorial is:

```bash
cd ~/nfds-tutorial/envs
```

To create the conda environments for the read cleaning step:

```bash
conda env create --file qc_assembly.yml
```

then activate the environment:

```bash
conda activate qc_assembly
```

Now you should be able to run the data pre-processing steps. For these examples we will be using sample [ERR065307](https://www.ncbi.nlm.nih.gov/sra/ERR065307) as a case study.

As previously mentioned, you can download both the raw read files and the adapter file needed to run trimmommatic (`NexteraPE-PE.fa`) from this [google drive folder](https://drive.google.com/drive/folders/1sSb7NM1VYhXeSLlXKdYsRk5M3yqsJlm0?usp=sharing).

These three files should be downloaded into the `nfds-tutorial/exercises/data_preprocessing/reads` folder within this repository.

### Read cleaning with [Trimmomatic](https://github.com/timflutre/trimmomatic)

First we navigate to the folder that contains the example raw reads:

```bash
cd  ~/nfds-tutorial/exercises/data_preprocessing/reads
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
cd  ~/nfds-tutorial/exercises/data_preprocessing/trimmed_reads
```
Now, we can run the FastQC command:

```bash
fastqc -o ~/nfds-tutorial/exercises/data_preprocessing/fast_qc ERR065307_R1_trimmed.fastq.gz ERR065307_R2_trimmed.fastq.gz
```
And to summarize all FastQC reports into one, we run MultiQC:

```bash
multiqc ~/nfds-tutorial/exercises/data_preprocessing/fast_qc
```
Now you can assess the read quality accross all samples.

### Assembly with [Unicycler](https://github.com/rrwick/Unicycler?tab=readme-ov-file#quick-usage)

Assemblies are one of the most computationally intensive steps in the pipeline. Thus, they are usually run in high performance computing clusters. It generally requires a fair bit of memory, and the program will fail if it runs out of memory. Thus, I recommend submitting this command as a job in your favorite high performance computer cluster. In the [scripts](/scripts) folder I have included a bash script named `run_unicycler_nfds_example.sh` that runs unicycler for the one sample we are using as an example. This script is formatted specifically for the [FasRC](https://www.rc.fas.harvard.edu/) cluster at Harvard University, so make sure to update the format/paths/etc. It is just meant to serve as an example.

The bash script mentioned above allocated 36 GB of memory across 16 cores. This might be overkill for a single sample, but it did run. Of course if you are running multiple samples make sure to allocate memory appropiately. Using that script, the assembly took 12 minutes and 32 seconds to run.

You can also run it from the command line directly, if you've allocated enough memory. Below is an example of how to do that.

Unicycler is a bit sensitive to package versions, so we create and activate a different environment from the `unicycler.yml` file:

```bash
conda env create --file unicycler.yml
conda activate unicycler
```

Now, from the folder containing the trimmed reads you can run Unicycler:

```bash
unicycler -1 ERR065307_R1_trimmed.fastq.gz -2 ERR065307_R2_trimmed.fastq.gz -o ~/nfds-tutorial/exercises/data_preprocessing/assemblies/ERR065307_assembly
```

Unicycler generates many output files, and there will be one folder per sample. In this case our folder is `~/ERR065307_assembly` because that is what we specified. The resulting assembly will be the file named `assembly.fasta` in this output folder.

### Assembly quality control with [Quast](https://github.com/ablab/quast)

Now the assembly quality can be assessed using Quast:

Due to package conflicts, quast is run from a separate conda environment. To create and use the quast environment:

```bash
conda deactivate # If there is another active environment
conda env create --file quast.yml
conda activate quast
```
Afterwards, quast can be run:

```bash
quast -o ~/nfds-tutorial/exercises/data_preprocessing/quast_qc/ERR065307_quast --threads 4 ~/nfds-tutorial/exercises/data_preprocessing/assemblies/ERR065307_assembly/assembly.fasta
```
The 'report.tsv' file in the generated folder contains a breakdown of various assembly quality metrics for your sample. This can be used to exclude assemblies based on quality. The exclusion criteria we use in our workflow are:

i) An N50 less than 15 kb; or
ii) ≥500 contigs, indicating the genome was too segmented; or 
iii) A genome length <1.9 Mb or > 2.4 Mb

## Pangenome analysis

We can now move on the the pangenome analysis exercise. By definition, pangenome analyses require more than one sample. Thus, for this part of the exercise we will be using 3 samples in the [Navajo population](https://www-ncbi-nlm-nih-gov.ezp-prod1.hul.harvard.edu/bioproject/PRJEB8327) as an example for the first couple of steps in this exercise. Running these analyses on a linux cluster is recommended, since some architectures are not compatible with the packages.

First download the 3 fasta files we will use from this [google drive folder](https://drive.google.com/drive/folders/1969ZeuSHZK8Xfx1wmXalc4ttI_JEMC44?usp=sharing). 

These three files should be downloaded into the `/nfds-tutorial/exercises/pangenome_analysis/assemblies` folder within this repository.

After data download, the first step is to navigate to the exercise subfolder:

```bash
cd ~/nfds-tutorial/exercises/pangenome_analysis
```
### Annotate assemblies with [Prokka](https://github.com/tseemann/prokka)

Now, we create and activate the conda environment containing prokka:

```bash
conda env create --file prokka_env.yml
conda activate prokka_env
```
Now we run a bash script that loops through every individual assembly file in a folder and runs prokka:

```bash
# Assign input arguments to variables
ASSEMBLIES_DIR="~/nfds-tutorial/exercises/pangenome_analysis/prokka_roary/assemblies"
OUTPUT_DIR="~/nfds-tutorial/exercises/pangenome_analysis/annotations"

# Create output directory if it does not exist
mkdir -p $OUTPUT_DIR

# Loop through each assembly file in the assemblies directory
for ASSEMBLY in ${ASSEMBLIES_DIR}/*.fasta; do
  # Get the base name of the file (without path and extension)
  BASENAME=$(basename ${ASSEMBLY} .fasta)

  # Run Prokka
  prokka --outdir ${OUTPUT_DIR}/${BASENAME} --prefix ${BASENAME} ${ASSEMBLY}
done
```
Prokka will create a folder per sample in the `annotations` folder. The annotation files will the the ones named `{sample_name}.gff` within each subfolder. 

This step can also take a bit of time (few minutes per sample), so you can also create a bash script and submit as a job to a high computing cluster, as we had previously discussed. 

### Get pangenome gene definitions with [Roary](https://github.com/sanger-pathogens/Roary)

As we are accustomed to, we create and activate the Roary environment:

```bash
conda env create --file roary_env.yml
conda activate roary_env
```
Now, let's copy our three .gff files into the `nfds-tutorial/exercises/pangenome_analysis/roary` folder:

```bash
find ~/nfds-tutorial/exercises/pangenome_analysis/annotations -type f -name "*.gff" -exec cp {} ~/nfds-tutorial/exercises/pangenome_analysis/roary \;
```

This command goes into each prokka subfolder, finds the .gff files and then copies them into the roary folder.

Afterwards, in the designated folder containing all annotation files, we can proceed to run Roary:

```bash
# Navigate to roary folder
cd ~/nfds-tutorial/exercises/pangenome_analysis/roary

# Run roary
roary -o {1} -i 90 -e -n -z -v *.gff
```

Roary assumes you are located in the folder that contains all the .gff files. It should not take too long to run on these three samples (), but for "real" runs with hundreds/thousands of samples, I highly recommend submitting as a cluster job allocating a fair bit of memory. I include an example of a bash script to do this in the [scripts](/scripts) folder named `torun_roary_cdc.sh`. Remember the script must be run from the folder containing the .gff files. In this script I allocate 1000 GB across 34 cores and for ~9000 samples it takes around 12 hours to run. 

Here Roary is run on the default parameters, except for the identity threshold which is set at 90% (instead of 98%). Intermediate files are kept, but you can decide the parameters that work best for your specific purpose. For more information on Roary's parameters and output files I highly recommend readind [Roary's documentation](https://sanger-pathogens.github.io/Roary/) which is very thorough. 

### Clean Roary gene definitions with [CLARC](https://github.com/IndraGonz/CLARC) (obtain final presence absence matrix for accessory genes)

CLARC is specifically meant to refine Clusters of Orthologous Genes (COG) assignmnents from a multi-population Roary run. So, our previous example Roary run is not compatible with the tool.

Thus, for the rest of this tutorial we will use the relevant Roary output files from a pangenome analysis run with >8000 _S. pneumoniae_ genomes encompassing more than 8 distinct carriage datasets. First step is installation.

#### Installing CLARC

Now, we have gotten to the exciting part of installing CLARC! As it is in the beta testing stage, we first build the given CLARC environment to install all dependencies and then we install the tool by copying the CLARC repository directly from GitHub.

For the installation, make sure you've installed git. If you have not, you can easily do so with:

```bash
conda install git
```
Afterwards you can proceed to install CLARC:

```bash
# Clone repository from GitHub (the repository will be cloned in your current working directory)
git clone https://github.com/IndraGonz/CLARC.git

# Navigate to envs folder
cd CLARC/envs

# Create environment from yml file
conda env create --file clarc_env.yml

# Activate environment
conda activate clarc_env

# Navigate to main CLARC folder
cd CLARC

# Install CLARC
python setup.py install
```

To check that clarc was installed correctly just write:

```bash
clarc --version
```

If it shows the version then it has been installed correctly!

#### Running CLARC

CLARC uses two output files from Roary: The 'gene_presence_absence.csv' which contains a presence absence matrix with all the COGs identified by Roary and 'pan_genome_reference.fa' which is a fasta file containing the representative sequences of these genes. The representative sequence is the longest instance of that COG across the samples. The .csv file can be big (>1 GB), as it contains ALL COG definitions and not only core/accessory.

Additionally, CLARC needs a text file with the names of the samples in your population of interest. The linkage constraints will be calculated based on this population. This can be all location from a particular geographic location, for example. This file should be named 'needed_sample_names.txt'. In this example, we will be running CLARC to refine the COG definitions in the Navajo population (937 samples).

You can download the CLARC input files from this [google drive folder](https://drive.google.com/drive/folders/1cfJtgLX9w2DGSb90K965lf_sxNn3k36f?usp=sharing).

These three inputs should be put in the `/nfds-tutorial/exercises/pangenome_analysis/clarc/data` subfolder. 

Afterwards, clarc can be run from the terminal (within the environment):

```bash
# Run CLARC
clarc --input_dir ~/nfds-tutorial/exercises/pangenome_analysis/clarc/data --output_dir ~/nfds-tutorial/exercises/pangenome_analysis/clarc/clarc_output
```

## Strain typing 

### Strain typing using [PopPUNK](https://github.com/bacpop/PopPUNK)

To be consistent with worlwide efforts to study pneumococcus, we type each sample for its global pneumococcal sequence cluster (GPSC). This is a PopPUNK typing developed by the [global pneumococcal sequencing project (GPS)] (https://www.pneumogen.net/gps/). 

To do this, we input an external clustering file into PopPUNK, which contains information from all the samples in the GPS database. This is done by following the [instructions for GPSC typing](https://www.pneumogen.net/gps/#/training#command-line) outlined by GPS.

**NOTE:** The version of PopPUNK supported by the GPS is not compatible with the Mac M1 architecture. Linux is recommended.

To run PopPUNK with external clustering, we need the database files provided by the GPS. These can be found in this [google drive folder](https://drive.google.com/drive/folders/1NvBlf_KX_H3ZuYwrNiGPoD5kBRHJwBkk?usp=sharing). Download both the `GPS_v8_external_clusters.csv` file and the `GPS_v8_ref` folder and place them in the `~/nfds-tutorial/exercises/poppunk/databases` folder within this repository. When running the poppunk command, you need to specify the location of these files.

PopPUNK uses the genome assemblies (.fasta files) to perform the GPSC assignment. We will use the assemblies found in the `nfds-tutorial/exercises/pangenome_analysis/assemblies` folder for this exercise. 

As input, you must provide a text file where the first column is your sample's name and the second is the path of the assembly. This is called a q file (I am not sure why). The needed input text file for our exercise (`nfds_example_qfile.txt`) is already provided in the `nfds-tutorial/exercises/poppunk/data` folder within this repository. If you have copied the repository to a specific location of your computer/cluster, **make sure to update the paths in the nfds_example_qfile.txt file**.

Generating the qfile can be annoying, so for your convenience I've also included here a bash [script](/scripts/make_poppunk_qfile.sh) that generates the qfile automatically by specifying the folder in which your fasta files are located. We don't need it for this example, however. 

After you've downloaded/updated all files in their specified locations, we can proceed to run PopPUNK!

First we create and activate the environment:

```bash
# Create poppunk environment
conda env create --file poppunk_env.yml

# Activate poppunk environment
conda activate poppunk_env
```

And now we can execute the poppunk assign command:

```bash
# Run PopPUNK
poppunk_assign --db ~/nfds-tutorial/exercises/poppunk/databases/GPS_v8_ref --external-clustering ~/nfds-tutorial/exercises/poppunk/databases/GPS_v8_external_clusters.csv --query ~/nfds-tutorial/exercises/poppunk/data/nfds_example_qfile.txt --output ~/nfds-tutorial/exercises/poppunk/data
```

The general GPSC assignments will be located in the 'poppunk_typing_external_clusters.csv' output file in the `/nfds-tutorial/exercises/poppunk/data` folder.

## Quadratic programming (QP) model 

### Predicting post-vaccine population structure with quadratic programming (QP) model

We've done it! We have gotten to the last step of the NFDS pipeline: Predicting post-vaccine population structure using pre-vaccine gene frequencies.

#### Preparing the data

In this example, we will be using the Southwest US (Navajo) dataset again. In this dataset, there was sampling before and after vaccination. This allows us to divide the samples by year of collection and we can calculate the gene and strain frequencies before and after vaccination.

But before we do that, let's download the data we will need for this part! What we will be using is the presence absence matrix that is output from CLARC (`clarc_condensed_presence_absence.csv`) and a metadata file that contains the PopPUNK GPSC results along with some other information about the samples, including time of collection (`MMN-meta-clean.xlsx`). Both of these files can be found in this [google drive folder](https://drive.google.com/drive/folders/1KSEkmj6v0Wmz4aJlz1_gs8vvf9MHygUq?usp=drive_link). This is the last google drive folder, I promise.

Download the files into the `nfds-tutorial/exercises/qp_prediction/data` folder within this repository. And now you are ready to get jupyter set up.

**Note:** I provided the 'CLARC presence absence' output here for ease, but if you ran the CLARC step in this tutorial you should have generated this exact same file. You can use this file instead of downloading it, which should be in this path: `nfds-tutorial/exercises/pangenome_analysis/clarc/clarc_output/clarc_results/clarc_condensed_presence_absence.csv`. 

#### Setting up the jupyter notebook

The quadratic programming model is run from a Jupyter Notebook. Instructions on how to install jupyter-lab can be found [here](https://jupyter.org/install).

Additionally, the QP solver itself is run through an R integration, so I've created an environment that includes the dependencies needed to run this integration. You can install it from the envs folder as per ush:

```bash
# Create r to python environment
conda env create --file rpy-env.yml
```

Then we activate:

```bash
# Activate r to python environment
conda activate rpy-env
```

**Within the rpy environment**, launch jupyter lab:

```bash
jupyter-lab
```
Now with the jupyter lab interface you can navigate to the `/exercises/qp_prediction` folder in this repository and open the file named `navajo_qp-prediction.ipynb`

#### Running the prediction within the jupyter notebook

Now you can navigate through the Jupyter notebook where you can run each step individually in subsequent cells. At the end there should be a plot where the predicted strain frequencies after vaccination are compared to the frequencies that were _actually_ observed after vaccination.

## Conclusion

This is a bare bones overview of the workflow used to study negative frequency dependent selection in _S. pneumoniae_. I hope you enjoyed running this as much as I enjoyed making it! :)

I would love to hear any problems, thoughs, suggestions, hopes, dreams, etc. Feel free to contact me with any of that at igonzalezojeda@g.harvard.edu


