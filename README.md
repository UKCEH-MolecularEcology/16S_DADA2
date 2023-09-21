# 16S_DADA2
# About
- This repository can be used for the analyses of 16S sequencing data via `DADA2` using a `snakemake` workflow

# Set up
```
git clone git@github.com:UKCEH-MolecularEcology/16S_DADA2.git
```

- You might want to adjust some settings in the `config/config.yaml`, e.g.
the name of the `snakemake` `conda` environment, paths and number of cores/threads.

### Setting up `samples.tsv` and `metadata.tsv`

#### Step1:
- There is no other option but to create a metadata file outside of this workflow and import
- NOTE: the first column of the metadata file should be named, `"Sample"`
- Import file and place in the `config` folder

#### Step2:
- Edit `line 15`, i.e. path to the fastq.gz files in the `notes/samples_tsv.sh` script
- From the `16S_DADA2` folder run the following:
```
bash notes/samples_tsv.sh
```

## Conda
[Conda user guide](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html)
```bash
# install miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod u+x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh # follow the instructions
# create venv
conda env create -f requirements.yml -n "YourEnvName"
```

# Run the analysis
```bash
conda activate snakemake
snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --cores 24 --conda-frontend conda --conda-prefix /hdd0/susbus/tools/conda_env --rerun-trigger mtime -rp
```
