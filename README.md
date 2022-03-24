# sbx_demic

[Sunbeam](https://github.com/sunbeam-labs/sunbeam) extension for [DEMIC](https://sourceforge.net/projects/demic/files/)

*Important* -- This was tested with modified demic scripts version 1.0.2, please check the sourceforge.net site for any updates. One change I made was to save a couple of R data files for "manual" analysis in Rstudio if one so desires. See the git history for any other changes.

## Installation

0. Install sbx_coassembly from https://github.com/sunbeam-labs/sbx_coassembly
1. git clone https://github.com/sunbeam-labs/sbx_demic
2. cp sbx_demic $SUNBEAM_DIR/extensions/
3. cat sunbeam/extensions/sbx_demic/config.yml >> sunbeam_config.yml (the config.yml that your are using for your given project)
4. Edit sunbeam_config.yml to have desired parameters
    - Make sure that
    - *all.paired_end* is true if you have paired end reads
    - *sbx_demic.extras* has any parameters you want to pass to DEMIC.pl

## Running

1. sunbeam run all_demic --use-conda --configfile */path/to/config* {rest of arguments for sunbeam} -- {arguments for snakemake}

### Trouble-shooting

- If you have trouble running with "--use-conda" it may be best just to install the needed packages into the sunbeam environment (e.g. `conda activate sunbeam && conda install --file sbx_demic_env.yml`) and then run sunbeam in two pieces (`sunbeam run all_coassemble --use-conda --configfile /path/to/config {rest of arguments for sunbeam}` and then `sunbeam run all_demic --configfile /path/to/config {rest of arguments for sunbeam}`). Please also create a new issue on GitHub detailing the error.

## References

- Sunbeam: https://github.com/sunbeam-labs/sunbeam

- DEMIC software: https://sourceforge.net/projects/demic/files/

- DEMIC publication: https://www.nature.com/articles/s41592-018-0182-0

- isolated Conda environment: http://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management
