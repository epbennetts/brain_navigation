#!/bin/csh
# Set name of job
#PBS -N DTBrainGenes
#PBS -o out_AccVsNumGenes_AllAreas.txt
#PBS -j oe
# Specify a queue:
# PBS -q physics
#PBS -l select=1:ncpus=12:mem=32GB
# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#PBS -l walltime=150:00:00
# Email user if job ends or aborts(and when it starts)
#PBS -m bea
#PBS -M eper8035@uni.sydney.edu.au
#PBS -V

# your commands/programs start here, for example:
cd "$PBS_O_WORKDIR"
# Show the host on which the job ran
hostname
# Launch the Matlab job
module load Matlab2019b
matlab -nodisplay -r "DT_AccVsNumGenes_ALL_AREAS; exit"
