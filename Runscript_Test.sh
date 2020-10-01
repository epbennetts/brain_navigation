#!/bin/csh
# Set name of job
#PBS -N DTBrainGenes
#PBS -o output.txt
#PBS -j oe
# Specify a queue:
# PBS -q physics
#PBS -l select=1:ncpus=1:mem=8GB
# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#PBS -l walltime=20:00:00
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
# matlab -nodisplay -singleCompThread -r "Test_Runscript; exit"
matlab -nodisplay -singleCompThread -r "Test; exit"
