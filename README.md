# nfds-tutorial
Simple tutorial of technical workflow to study negative frequency dependent selection (NFDS) in the accessory genome of Streptococcus pneumoniae. Covers the bioinformatic steps to go from shotgun whole genome sequencing raw reads to accessory gene frequencies. The accessory gene frequencies can be used to predict the post-vaccination population structure of Streptococcus pneumoniae. The implementation of this NFDS-based quadratic programming predictive model is also discussed.

# Step 1: Trim reads and assemble with Unicycler

# Step 2: Assembly QC with Quast

# Step 3: Obtain accessory gene definitions

## Step 3a: Annotate assemblies with prokka

## Step 3b: Get pangenome gene definitions with Roary

## Step 3c: Clean Roary gene definitions with CLARC (obtain final presence absence matrix for accessory genes)

# Step 4: Strain typing using poppunk

# Step 5: Predicting post-vaccine population structure with quadratic programming (QP) model


