# sbx_vsearch

[Sunbeam] extension for [DEMIC] (https://sourceforge.net/projects/demic/files/)

*Important* The demic scripts under /demic are v1.0.2, please check the sourceforge.net site for any updates 

## Installation

1. git clone https://github.com/sunbeam-labs/sbx_demic
2. cp sbx_demic $SUNBEAM_DIR/extensions/
3. cat sunbeam/extensions/sbx_demic/config.yml >> sunbeam_config.yml (the config.yml that your are using for your given project)
4. edit sunbeam_config.yml to have desired parameters (important: contigs_dir points to your pre-assembled metagenome)
5. copy demic scripts ([DEMIC.pl](/demic/DEMIC.pl), estGrowthRate1.0.2.R, and testR1.0.2.R) to your local bin directory (e.g. $HOME/bin)

## Running

1. sunbeam all_demic {rest of arguments for sunbeam}

## References

[Sunbeam]: https://github.com/sunbeam-labs/sunbeam
[DEMIC software]: https://sourceforge.net/projects/demic/files/
[DEMIC publication]: https://www.nature.com/articles/s41592-018-0182-0
[isolated Conda environment]: http://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management
