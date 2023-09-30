<img src="https://github.com/sunbeam-labs/sunbeam/blob/stable/docs/images/sunbeam_logo.gif" width=120, height=120 align="left" />

# sbx_demic

<!-- Begin badges -->
[![Tests](https://github.com/sunbeam-labs/sbx_demic/actions/workflows/test.yml/badge.svg)](https://github.com/sunbeam-labs/sbx_demic/actions/workflows/test.yml)
[![Super-Linter](https://github.com/sunbeam-labs/sbx_demic/actions/workflows/linter.yml/badge.svg)](https://github.com/sunbeam-labs/sbx_demic/actions/workflows/linter.yml)
[![Conda Envs Status](https://byob.yarr.is/sunbeam-labs/sbx_demic/env_check)](https://github.com/sunbeam-labs/sbx_demic/actions/workflows/check_conda_envs.yml)
[![DOI:10.1038/s41592-018-0182-0](https://badgen.net/badge/Published%20in/Nat%20Methods/blue)](https://doi.org/10.1038/s41592-018-0182-0)
<!-- End badges -->

A [Sunbeam](https://github.com/sunbeam-labs/sunbeam) extension for estimating bacterial growth rates via peak-to-trough ratios (PTRs) using [DEMIC](https://sourceforge.net/projects/demic/files/). In preparation to use demic, reads are first assembled using [Megahit](https://github.com/voutcn/megahit), then binned by inferred genome using [MaxBin2](https://sourceforge.net/projects/maxbin2/), after which reads are mapped back onto contigs using [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and [Samtools](https://github.com/samtools/samtools).

Note: The original paper for demic can be found [here](https://doi.org/10.1038/s41592-018-0182-0) and the associated software [here](https://sourceforge.net/projects/demic/files/). A lot of work has gone into demic since publication to make it more reliable, installable, and easy to use. The main functionality of demic is now housed in a python library for generating special coverage files ([pycov3](https://github.com/Ulthran/pycov3)) and an R package for estimating PTRs from these coverage files ([DEMIC](https://github.com/Ulthran/DEMIC/tree/master)).

## Installation

To install, activate your conda environment (using the name of your environment) and use `sunbeam extend`:

    conda activate <i>sunbeamX.X.X</i>
    sunbeam extend https://github.com/sunbeam-labs/sbx_demic.git

## Usage

To generate PTRs, create a project and use the `all_demic` target:

    sunbeam init --data_fp /path/to/reads/ /path/to/project/
    sunbeam run --profile /path/to/project/ all_demic

## Configuration

  - coassembly_threads: Is the number of threads to use in assembly steps
  - group_file: Is the path to the coassembly group file, leave blank to treat all samples as one group
  - demic_threads: Is the number of threads to use in the binning, mapping, and estimation steps
  - single_genome: Should be true if your samples only contain a single genome, false otherwise
  - extras: Are parameters to be passed to PyCov3

Note: See [sbx_coassembly](https://github.com/sunbeam-labs/sbx_coassembly) README for more information on how to group samples with the group_file. See [PyCov3](https://github.com/Ulthran/pycov3) for more information on what parameters can be passed for coverage generation.

## Legacy Installation

For sunbeam versions <3 or if `sunbeam extend` isn't working, you can use `git` directly to install an extension:

    git clone https://github.com/sunbeam-labs/sbx_demic.git extensions/sbx_demic

and then include it in the config for any given project with:

    cat extensions/sbx_demic/config.yml >> /path/to/project/sunbeam_config.yml
