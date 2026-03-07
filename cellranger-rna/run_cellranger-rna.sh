#!/bin/bash

#SBATCH --ntasks=1 --cpus-per-task=64 
#SBATCH --time=50:00:00
#SBATCH --job-name=cellranger-rna
***REMOVED***

#this is a simple bash script to run the cellranger-rna in the Vlaams supercomputer.

#parameters
INPUTS_DIR="/scratch/antwerpen/205/vsc20542/cellranger-rna/input_fastq/"
OUTPUTS_DIR="/scratch/antwerpen/205/vsc20542/cellranger-rna/cellranger_outputs/ran_with_LdCL_genome/"
CELLRANGER_DIR="/scratch/antwerpen/205/vsc20542/tools/cellranger-rna/cellranger-8.0.1/"
TRANSCR_REF="/scratch/antwerpen/grp/aitg/jcdujardin/reference_genomes/L_donovani/LdCL/cellrangerRNA_formatted/LdCL_no_multicopy"
SAMPLE="105875-001-002"
CREATE_BAM=true

#attempt to remove output folder for the sample (in case previous run was made)
rm -rf "$OUTPUTS_DIR/$SAMPLE"

#attempt to create general output folder in case it doesn't exist yet. 
mkdir -p "$OUTPUTS_DIR"

#get the modules
module --force purge
module load calcua/all

# run cellranger
$CELLRANGER_DIR/cellranger count \
    --id="$SAMPLE" \
    --output-dir="$OUTPUTS_DIR/$SAMPLE" \
    --transcriptome="$TRANSCR_REF" \
    --fastqs="$INPUTS_DIR" \
    --sample="$SAMPLE" \
    --create-bam=$CREATE_BAM \
    --r1-length=28 \
    --r2-length=91

# after cellranger is ran, move the outputs one level up
mv $OUTPUTS_DIR/$SAMPLE/outs/* $OUTPUTS_DIR$SAMPLE
rm -r $OUTPUTS_DIR/$SAMPLE/outs
