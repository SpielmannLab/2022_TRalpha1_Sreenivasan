#! /bin/bash
# Run this script to invoke multiple sbatch for parallel processing of cellranger count

# First declare directories where all the data is
path_fastqdata="/data/humangen_mouse/hypthmsMittag/data"
path_refdata="/data/humangen_mouse/refdata/refdata-gex-mm10-2020-A"
path_results="/data/humangen_mouse/hypthmsMittag/results"
echo "Message: All paths have been set"

# Define samples in a array
samples=( "mpimg_L20986-1_B-Tra1" "mpimg_L20987-1_C-wt" )
samples=( "mpimg_L23584-1_Tra1-mutant" "mpimg_L23585-1_Tra1-WT")


for sample in ${samples[@]}
do
	sbatch --job-name=count_${sample} script_count_core.sh $sample $path_refdata $path_fastqdata $path_results
	echo "Message: Job count_${sample} has been sent to script_count_core.sh"	
done

echo "Message: All samples have been sent for count pipeline"
