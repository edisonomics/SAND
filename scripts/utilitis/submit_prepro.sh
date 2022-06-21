#! /bin/bash
#SBATCH --job-name=fid_prepro        # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=120gb                     # Job memory request
#SBATCH --time=150:00:00               # Time limit hrs:min:sec
#SBATCH --output=fid_prepro.%j.out     # Standard output log
#SBATCH --error=fid_prepro.%j.err      # Standard error log
#SBATCH --cpus-per-task=45             # Number of CPU cores per task
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Your@email  # Where to send mail

cd $SLURM_SUBMIT_DIR

module load  matlab/R2020b

#local folder for matlab temp files
mkdir temp_prepro
mkdir -p ./temp_prepro/$SLURM_JOB_ID
time matlab -nodisplay < data_preprocess.m > matlab_${PBS_JOBID}.out
