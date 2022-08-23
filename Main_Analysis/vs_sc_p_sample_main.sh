#! /bin/bash
# Run this script to carry out the entire Single-Cell-RNA analysis pipeline
# STEP1. Make sure you have run the mkfastq_main.sh and count_main.sh scripts
# STEP2. Provide all the details in the file "vs_sc_p_sample.yaml"
# STEP3. Close this file and execute it using ./vs_sc_p_sample_main.sh

sbatch --job-name=vs_sc_p_sample /data/humangen_mouse/scpipeline/src/vs_sc_p_sample_core.sh 
