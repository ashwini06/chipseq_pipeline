#! /bin/bash -l
#SBATCH -A b2012025
#SBATCH -M milou
#SBATCH -p core -n 1
#SBATCH -t 1:00:00
#SBATCH -J chipPipeline
#SBATCH -e /proj/b2012025/RAW_DATA/test/mapping.stderr.txt
#SBATCH -o /proj/b2012025/RAW_DATA/test/mapping.stdout.txt
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ashwini.jeggari@scilifelab.se


cd /proj/b2012025/RAW_DATA/test
# ==============================================================================
# LOADING MODULES
module load bioinfo-tools
source activate master
cd /proj/b2012025/RAW_DATA/test

# RUNNING SNAKEMAKE MAPPING SCRIPT WITH 16 cores
#snakemake --snakefile Snakefile 
snakemake -R find_narrow_peaks --unlock 
