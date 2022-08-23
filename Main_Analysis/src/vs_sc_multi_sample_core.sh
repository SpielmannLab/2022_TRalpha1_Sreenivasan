#! /bin/bash

### Submit this Script with: sbatch <script.sh> ###

# Parameters for slurm (don't remove the # in front of #SBATCH!)
#  Use partition shortterm/debug/longterm:
#SBATCH --partition=shortterm
#  Use so many node:
#SBATCH --nodes=1
#  Request so many cores (hard constraint):
#SBATCH -c 8
#  Request so much of memory (hard constraint):
#SBATCH --mem=400GB
# set slurm file output nomenclature
#SBATCH --output "slurm-%x-%j.out"

PATH=$WORK/.omics/anaconda3/bin:$PATH #add the anaconda installation path to the bash path
source $WORK/.omics/anaconda3/etc/profile.d/conda.sh # some reason conda commands are not added by default

# Load your necessary modules:
conda activate ScRNA

# Move to SCRATCH were everything will be done
cd $SCRATCH

# Submit the Nextflow Script:
nextflow run /data/humangen_mouse/scpipeline/src/sc_multi_sample.nf -params-file /data/humangen_mouse/scpipeline/Varun/vs_sc_multi_sample.yaml --id ${SCRATCH/"/scratch/"/}
