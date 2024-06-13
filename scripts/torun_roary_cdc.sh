#! /bin/bash
#SBATCH -p hsph
#SBATCH -N 1
#SBATCH -n 34
#SBATCH -t 3-00:00
#SBATCH --mem=1000000
#SBATCH -o roary.%N.%j.out
#SBATCH -e roary.%N.%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=igonzalezojeda@g.harvard.edu

export OMP_NUM_THREADS=$SLURM_NTASKS

# Activate roary environments
source activate roary_env

#Run roary for pangenome analysis
roary -p $SLURM_NTASKS -o {1} -i 90 -e -n -z -v *.gff
