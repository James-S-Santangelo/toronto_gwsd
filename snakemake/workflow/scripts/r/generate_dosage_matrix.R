library(algatr)
library(readr)
library(tibble)

vcf <- vcfR::read.vcfR(snakemake@input[["vcf"]])
dos_mat <- as.data.frame(t(vcf_to_dosage(vcf))) %>%
    rownames_to_column("site")
write_delim(dos_mat, snakemake@output[[1]], delim="\t")
