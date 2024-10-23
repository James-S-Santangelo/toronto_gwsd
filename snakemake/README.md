## Analyses of genomic data for Toronto_GWSD

### Description of repository

This repository contains code necessary to reproduce the genomic analyses in
the Toronto_GWSD manuscript. Raw reads from this project have been deposited in
GenBank (BioProject #XXX). This pipeline runs from the raw reads and includes
read trimming, mapping, and QC of reads and alignments, as well as variant
calling, population structure analyses, and analyses of selective sweeps using
SFS and haplotype-based approaches. The pipeline uses `Conda`, `Snakemake`, and
`Singularity` for workflow management and reproducibility. All Snakefiles and
directories are well-documented, but here is a brief overview of the pipeline
and directories in this repository:

#### Overview of directories

- [config](./config): Snakemake configuration files for different clusters.
- [notebooks](./notebooks): Jupyter Notebooks detailing analyses of diversity,
  population structure, and HCN/Ac/Li differentiation.
- [resources](./resources): Text files used in pipeline (e.g., sample
  information, chromosomes, etc.)
- [workflow](./workflow): Main Snakemake workflow with rules, environments,
  scripts, notebooks, and cluster profiles for running the pipeline on Compute
  Canada SLURM-based clusters.

### Using the pipeline

This pipeline requires `Conda` and `Singularity`:

- A minimal installation of `Conda` (i.e., Miniconda) can be installed by
  following the instructions for your platform
  [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
- Installation of `Singularity` requires Admin privileges, but using
  `Singularity` to run pre-created containers does not. Installation
  instructions can be found
  [here](https://sylabs.io/guides/3.5/admin-guide/installation.html). All
  `Singularity` containers used in this pipeline are avalaible in [this public
  reposity](https://cloud.sylabs.io/library/james-s-santangelo), though they
  will be automatically pulled and executed by the pipeline. 

Assuming `Conda` is installed, the this repository's `Conda` environment can be
replicated by running the following command:

`conda env create -f environment.yaml -n gwsd`

This will create a `Conda` environment named _gwsd_ containing a minimal
set of dependencies required to run the pipeline (e.g., Python 3 and
Snakemake 7.21).

After activating the environment (`conda activate gwsd`), the pipeline can
be executed from the [workflow](./workflow) directory by running a command that
looks something like:

`snakemake --use-conda --use-singularity --singularity-args "--bind <path>"
--configfile ../config/<configfile> --notemp -j <cores>`

for local execution. Here, `<path>` is the path on the cluster from which files
will be read/written (e.g., `/scratch`), `<configfile>` is one of the
configfiles in the [config](./config) directory that needs to be modified
to match the paths on your system, and `<cores>` is the number of cores
available for executing parallel processes. 

For execution on a SLURM cluster, the pipeline can be executed by running:

`snakemake --profile compute-canada --configfile ../config/<configfile>`

Note that the YAML configfiles in the
[compute-canada](./workflow/compute-canada/) directory will likely need to be
modified to accomodate the paths on your cluster.
