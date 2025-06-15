[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14014248.svg)](https://doi.org/10.5281/zenodo.14014248)

## Genome-Wide Selection and Demography in Toronto white clover populations (Toronto_GWSD) 

### Manuscript: Signatures of selective sweeps in urban and rural white clover populations

### Abstract

Urbanization is increasingly recognized as a powerful force of evolutionary
change. However, anthropogenic sources of selection can often be similarly
strong and multifarious in rural habitats, and whether selection differs in
either strength or its targets between habitats is rarely considered. Despite
numerous examples of phenotypic differentiation between urban and rural
populations, we still lack an understanding of the genes enabling adaptation to
these contrasting habitats. In this study, we conducted whole genome sequencing
of 120 urban, suburban, and rural white clover plants from Toronto, Canada and
used these data to identify urban and rural signatures of positive selection. We
found evidence for selection in genomic regions involved in abiotic stress
tolerance and growth/development in both urban and rural populations, and clinal
change in allele frequencies at SNPs within these regions. Patterns of allele
frequency and haplotype differentiation suggest that most sweeps are incomplete
and our strongest signals of selective sweeps overlap known large-effect
structural variants. These results highlight how both urban and rural habitats
are driving ongoing selection in Toronto white clover populations, and motivate
future work disentangling the genetic architecture of ecologically important
phenotypes underlying adaptation to contemporary anthropogenic habitats.

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
