import tskit
import pandas as pd

results = {"num_sites": [], "num_trees": []}
region = snakemake.wildcards["n"]
for s in range(snakemake.params["n_samples"]):
    ts = tskit.load(f"{snakemake.params['arg_path']}/trees/region{region}/region{region}_{s}.trees")
    results["num_sites"].append(ts.num_sites)   
    results["num_trees"].append(ts.num_trees) 
df_out = pd.DataFrame(results)
df_out["regionID"] = region
df_out.to_csv(snakemake.output[0], index = False)
