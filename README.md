[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14014248.svg)](https://doi.org/10.5281/zenodo.14014248)

## Genome-Wide Selection and Demography in Toronto white clover populations (Toronto_GWSD) 

### Manuscript: Signatures of selective sweeps in urban and rural white clover populations

### Abstract

Urbanization is increasingly recognized as a powerful force of evolutionary
change. However, anthropogenic sources of selection can often be similarly
strong and multifarious in rural habitats, but these are often ignored in
studies of urban evolutionary ecology. Despite numerous examples of phenotypic
differentiation between urban and rural populations, we still lack an
understanding of the genes enabling adaptation to these contrasting habitats
and the  genetic architecture underlying urban and rural adaptation. In this
study, we conducted whole genome sequencing of 120 urban, suburban, and rural
white clover plants from Toronto, Canada. We used these data to identify
signatures of selective sweeps across the genome using both SFS and
haplotype-based approaches, and characterize the architecture of selective
sweeps. We found evidence for selection in genomic regions involved in abiotic
stress tolerance and growth/development in both urban and rural populations.
Urban and rural populations did not differ in the proportion of hard vs. soft
sweeps, though urban populations were characterized by wider sweeps, which may
indicate differences in the strength or timescale of selection. In addition,
patterns of allele frequency and haplotype differentiation suggest that most
sweeps are incomplete. These results highlight how both urban and rural
habitats are driving ongoing selection in white clover populations, and
motivate future work disentangling the genetic architecture of ecologically
important phenotypes, and estimating the strength and timescale of selection
underlying adaptation to contemporary anthropogenic habitats.

### Description of repository

This repository contains code and data necessary to reproduce the manuscript's
results. The repo can be clones using the following command:

`git clone https://github.com/James-S-Santangelo/toronto_gwsd.git`

Here is a brief description of each subdirectory. Details and documentation can
be found in each subdirectory:

- [figures_tables](./figures_tables): Contains figures and tables generated as
  part of the Snakemake pipeline below
- [snakemake](./snakemake): Contains the pipeline used to generate the genomic
  results Uses Conda + Snakemake for reproducibility and pipeline management. 
