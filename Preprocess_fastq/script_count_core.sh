#! /bin/bash
### Submit this Script with: sbatch script.sh ###
 
# Parameters for slurm (don't remove the # in front of #SBATCH!)
#  Use partition debug:
#SBATCH --partition=shortterm
#  Use one node:
#SBATCH --nodes=1
#  Request 32 cores (hard constraint):
#SBATCH -c 32
#  Request 128GB of memory (hard constraint):
#SBATCH --mem=128GB
#  Notify me at job start and end:
#SBATCH --mail-type=ALL
#  Send the notifications to (change to your email address!):
#SBATCH --mail-user=varun.sreenivasan@uni-luebeck.de
#  Find your job easier with a name:
#SBATCH --job-name=cellranger_count

if [ $# -eq 4 ]; then
	sample="$1"
	path_refdata="$2"
	path_fastqdata="$3"
	path_results="$4"
	echo "Message: All Arguments have been Read for the job"
else
	echo "Usage: $0 <sample> <path_refdata> <path_fastqdata> <path_results>"
	exit
fi

# Load your necessary modules (example):
module load cellranger/5.0.1

#Prepare $SCRATCH
mkdir -p  $SCRATCH/inputdata/{fastq,refdata} $SCRATCH/scripts/ $SCRATCH/results/

# Copy input data into your local job directory for much faster execution:
cp -v $path_fastqdata/${sample}*.gz $SCRATCH/inputdata/fastq/
cp -av ${path_refdata}/* $SCRATCH/inputdata/refdata/
echo "Data copied to scratch successfully"

# Set path to the results, where cell ranger will output the data
cd "$SCRATCH/results/"
echo "The current directory has been set to $(pwd)"

# Run cellranger count routine
# cellranger_5.0.1.simg cellranger count --id=$sample \
# 	--fastqs="$SCRATCH/inputdata/fastq" \
# 	--transcriptome="$SCRATCH/inputdata/refdata" \
# 	--sample=$sample

# Run cellranger count including introns routine
cellranger_5.0.1.simg cellranger count --include-introns \
	--id="${sample}_wIntrons" \
 	--fastqs="$SCRATCH/inputdata/fastq" \
 	--transcriptome="$SCRATCH/inputdata/refdata" \
 	--sample=$sample

echo "My message: Cell Ranger execution completed"

# Save the results from the scratch folder
cp -a $SCRATCH/results/* $path_results/
