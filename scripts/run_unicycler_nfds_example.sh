#! /bin/bash
#SBATCH -p hsph
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 1-00:00
#SBATCH --mem=64000
#SBATCH -o unicycler.%N.%j.out
#SBATCH -e unicycler.%N.%j.err
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=igonzalezojeda@g.harvard.edu

export OMP_NUM_THREADS=$SLURM_NTASKS

# Activate environment
source activate unicycler

#Run unicycler for assembly
unicycler -1 /n/holylfs05/LABS/lipsitch_lab/Lab/hsphs10/igonzalezojeda/strep_nfds/CDC_NFDS_training/exercises/data_preprocessing/trimmed_reads/ERR065307_R1_trimmed.fastq.gz -2 /n/holylfs05/LABS/lipsitch_lab/Lab/hsphs10/igonzalezojeda/strep_nfds/CDC_NFDS_training/exercises/data_preprocessing/trimmed_reads/ERR065307_R2_trimmed.fastq.gz -o /n/holylfs05/LABS/lipsitch_lab/Lab/hsphs10/igonzalezojeda/strep_nfds/CDC_NFDS_training/exercises/data_preprocessing/assemblies/ERR065307_assembly
