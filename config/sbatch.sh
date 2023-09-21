#!/bin/bash -l

# slurm settings if called using sbatch
#SBATCH -J CIRCLES
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 15
#SBATCH --time=0-06:00:00
#SBATCH -p batch
#SBATCH --qos=normal

# SMK_JOBS=${SLURM_CPUS_PER_TASK:-1} # get from slurm or set to 1
SMK_JOBS=48
SMK_ENV="snakemake" # snakemake conda env (see requirements.yaml)
SMK_PRF="/mnt/data/sbusi/vip/conda" # snakemake conda prefix
SMK_CONFIG="config/config.yaml" # sankemake config file

# Checks
echo "conda activate: ${SMK_ENV}"
conda activate ${SMK_ENV}
if [[ -z "${SMK_PRF}" ]]; then 
    echo "Error: No snakemake conda prefix given"
    exit 1
elif [[ ! -d "${SMK_PRF}" ]]; then 
    echo "Error: ${SMK_PRF} does not exist"
    exit 1
else 
    echo "snakemake conda prefix: ${SMK_PRF}"
fi

# snakemake -s workflow/Snakefile -rpn --jobs 1 --local-cores 1 --unlock --configfile "${SMK_CONFIG}"
snakemake -s workflow/Snakefile -rp --jobs ${SMK_JOBS} --local-cores 1 --rerun-incomplete \
--configfile "${SMK_CONFIG}" --use-conda --conda-prefix ${SMK_PRF}
