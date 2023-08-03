#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=200GB
#SBATCH --cpus-per-task 3
#SBATCH --time=120:00:00
#SBATCH --output=/data/davis_lab/sealockj/projects/mdd_wbc_long/ad_response/add_antipsychotics/sd_pull_20230511/ad_response/pull_ad_out.txt
#SBATCH --mail-user=jsealock@broadinstitute.org
#SBATCH --mail-type=ALL
#SBATCH --job-name="pull_ad"
#SBATCH --account=davis_lab

module load GCC/8.2.0  OpenMPI/3.1.4 R/3.6.0
cd /data/davis_lab/sealockj/projects/mdd_wbc_long/ad_response/add_antipsychotics/sd_pull_20230511/ad_response/

Rscript get_antidepressants_0511_pull.R 

nohup Rscript get_antidepressants_20230607_pull.R &


