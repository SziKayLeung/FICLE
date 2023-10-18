#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=1:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address


##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2
source activate ficle
FICLE_ROOT=/lustre/projects/Research_Project-MRC148213/sl693/scripts/FICLE
export PATH=$PATH:${FICLE_ROOT}


##-------------------------------------------------------------------------

# directory paths
inputDir=/lustre/projects/Research_Project-MRC148213/sl693/scripts/FICLE/vignette/0_input/
outputDir=/lustre/projects/Research_Project-MRC148213/sl693/scripts/FICLE/vignette/0_results/

# run ficle
ficle.py --gene=Trem2 \
    --reference=${inputDir}/subsetted_gencode_vM22_annotation.gtf \
    --input_bed=${inputDir}/input_sorted.bed12 \
    --input_gtf=${inputDir}/input.gtf  \
    --input_class=${inputDir}/input_classification.txt  \
    --output_dir=${outputDir}
    