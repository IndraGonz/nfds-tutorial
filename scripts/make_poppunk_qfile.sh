#!/bin/bash

# Define the directory containing the fasta files
FASTA_DIR="nfds-tutorial/exercises/pangenome_analysis/assemblies"

# Define the output file
OUTPUT_FILE="nfds-tutorial/exercises/poppunk/data/nfds_example_qfile.txt"

# Clear the output file if it already exists
> "$OUTPUT_FILE"

# Loop through each fasta file in the directory
for FILE in "$FASTA_DIR"/*.fasta; do
  # Get the file name without the extension
  BASENAME=$(basename "$FILE" .fasta)
  # Write the file name and full path to the output file
  echo -e "${BASENAME}\t${FILE}" >> "$OUTPUT_FILE"
done
