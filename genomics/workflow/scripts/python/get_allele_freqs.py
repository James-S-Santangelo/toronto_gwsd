import cyvcf2
import pandas as pd

sample_df = pd.read_table(snakemake.input["samples"])
sample_df = sample_df[["Sample", "Habitat"]]

bams = pd.read_table(snakemake.input["bams"][0], names=["bam"])
bams["Sample"] = bams["bam"].str.extract("(s_\d+_\d+)")
bams = bams.merge(sample_df, on = "Sample", how = "left")[["Sample", "Habitat"]]

urban_samples = bams[bams["Habitat"] == "Urban"]["Sample"].tolist()
rural_samples = bams[bams["Habitat"] == "Rural"]["Sample"].tolist()

urban_vcf = cyvcf2.VCF(snakemake.input["vcf"], samples=urban_samples)
rural_vcf = cyvcf2.VCF(snakemake.input["vcf"], samples=rural_samples)

af_dict = {"chrom" : [],
           "pos" : [],
           "rural_af" : [],
           "urban_af" : []}

for variant in urban_vcf:
    af_dict["chrom"].append(variant.CHROM)
    af_dict["pos"].append(variant.POS)
    af_dict["urban_af"].append(variant.aaf)

for i, variant in enumerate(rural_vcf):
    assert af_dict["chrom"][i] == variant.CHROM
    assert af_dict["pos"][i] == variant.POS
    af_dict["rural_af"].append(variant.aaf)

df = pd.DataFrame(af_dict)
df.to_csv(snakemake.output[0], sep="\t", index=False)
