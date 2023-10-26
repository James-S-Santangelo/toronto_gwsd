import tskit
import pandas as pd

# Get sample order in VCF and habitats
sample_df = pd.read_csv(snakemake.config["samples"], sep='\t')
sample_df = sample_df[["Sample", "Habitat"]]
bams = pd.read_table(snakemake.input["bams"][0], names = ["bam"])
bams["Sample"] = bams["bam"].str.extract("(s_\\d+_\\d+)")
bams = bams.merge(sample_df, on="Sample", how="left")[["Sample", "Habitat"]]
bams_filt = bams[bams["Habitat"].isin(["Urban", "Rural"])].reset_index()

# Get sample and haplotype indices
urban_sample_indices = bams_filt.index[bams_filt['Habitat'] == "Urban"].tolist()
rural_sample_indices = bams_filt.index[bams_filt['Habitat'] == "Rural"].tolist()
urban_hap_indices = []
for i in urban_sample_indices:
    hap1_i = i*2
    hap2_i = (i*2) + 1
    urban_hap_indices.append(hap1_i)
    urban_hap_indices.append(hap2_i)

rural_hap_indices = []
for i in rural_sample_indices:
    hap1_i = i*2
    hap2_i = (i*2) + 1
    rural_hap_indices.append(hap1_i)
    rural_hap_indices.append(hap2_i)

# Load trees using tskit and calculate site and branch-based Fst
iter_list = []
site_fsts = []
branch_fsts = []
for tree in snakemake.input["trees"]:
    region = tree.split("/")[8]
    iter = tree.split("/")[9].split("_")[1].split(".")[0]
    iter_list.append(iter)
    ts = tskit.load(tree)
    urban_haps = ts.samples()[urban_hap_indices]
    rural_haps = ts.samples()[rural_hap_indices]
    site_fst = ts.Fst([urban_haps, rural_haps], mode="site")
    site_fsts.append(site_fst)
    branch_fst = ts.Fst([urban_haps, rural_haps], mode="branch")
    branch_fsts.append(branch_fst)

# Convert Fsts to dataframe and write to disk
fst_df = pd.DataFrame(zip(site_fsts, branch_fsts, iter_list), columns=["site", "branch", "iter"])
fst_df_mod = fst_df.melt(id_vars="iter", var_name="type", value_name="fst")
fst_df_mod["regionID"] = snakemake.wildcards["n"]
fst_df_mod.to_csv(snakemake.output[0], sep="\t", index=False)
