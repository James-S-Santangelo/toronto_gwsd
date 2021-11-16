# Python functions used throughout snakemake workflow

def get_raw_reads(wildcards):
    """
    Extract forward and reverse read FASTQ paths from file
    """
    R1 = glob.glob(RAW_READ_DIR + '/{0}/{0}_*_1.fq.gz'.format(wildcards.sample))[0]
    R2 = glob.glob(RAW_READ_DIR + '/{0}/{0}_*_2.fq.gz'.format(wildcards.sample))[0]
    return { 'read1' : R1, 'read2' : R2 }

def get_fastas_to_concat(wildcards):
    if wildcards.gene == 'rbcl':
        return expand(rules.chloroplast_gene_consensus.output, sample=SAMPLES, gene='rbcl')
    elif wildcards.gene == 'matk':
        return expand(rules.chloroplast_gene_consensus.output, sample=SAMPLES, gene='matk')

def aggregate_ncbi_input(wildcards):
    checkpoint_output = checkpoints.download_nt_database.get(**wildcards).output[0]
    DB = glob_wildcards(os.path.join(checkpoint_output, 'nt.{db}.tar.gz')).db
    return expand('{0}/ncbi_nt_database/nt.{{db}}.{{ext}}'.format(PROGRAM_RESOURCE_DIR), db=DB, ext=['nhd', 'nhi', 'nhr', 'nin', 'nnd', 'nni', 'nog', 'nsq'])

def get_vcfs_by_chrom(wildcards):
    return expand(rules.freebayes_call_variants.output, chrom=wildcards.chrom, i=FREEBAYES_CHUNKS)

def get_bed_to_subset(wildcards):
    all_bed_files = rules.get_fourfold_zerofold.output
    bed = [bed for bed in all_bed_files if wildcards.site in os.path.basename(bed)]
    return bed

def angsd_sfs_input(wildcards):
    saf_idx = rules.angsd_saf_likelihood_allSites.output.saf_idx
    sites_idx = rules.angsd_index_sites.output.idx
    if wildcards.site == 'allSites':
        sites = rules.extract_angsd_allSites.output
    else:
        sites = rules.split_angsd_sites_byChrom.output
    return { 'saf_idx' : saf_idx, 'sites_idx' : sites_idx, 'sites' : sites }

def angsd_estimate_thetas_input(wildcards):
    saf_idx = rules.angsd_saf_likelihood_allSites.output.saf_idx
    sfs = rules.angsd_estimate_sfs.output
    sites_idx = rules.angsd_index_sites.output.idx
    if wildcards.site == 'allSites':
        sites = rules.extract_angsd_allSites.output
    else:
        sites = rules.split_angsd_sites_byChrom.output
    return { 'saf_idx' : saf_idx, 'sfs' : sfs, 'sites_idx' : sites_idx, 'sites' : sites }

def get_angsd_stats_toConcat(wildcards):
    return expand(rules.angsd_diversity_neutrality_stats.output, chrom=CHROMOSOMES, site=wildcards.site)

def get_angsd_sfs_toConcat(wildcards):
    return expand(rules.angsd_estimate_sfs.output, chrom=CHROMOSOMES, site=wildcards.site)

def get_sites_for_angsd_index(wildcards):
    if wildcards.site == 'allSites':
        return rules.extract_angsd_allSites.output
    elif wildcards.site == '0fold' or wildcards.site == '4fold':
        return rules.split_angsd_sites_byChrom.output

def get_angsd_gl_toConcat(wildcards):
    if wildcards.site == 'allSites':
        return expand(rules.angsd_gl_allSites.output.gls, chrom=CHROMOSOMES, maf=wildcards.maf, site=wildcards.site)
    else:
        return expand(rules.subset_angsd_gl.output, chrom=CHROMOSOMES, maf=wildcards.maf, site=wildcards.site)

def get_angsd_maf_toConcat(wildcards):
    if wildcards.site == 'allSites':
        return expand(rules.angsd_gl_allSites.output.mafs, chrom=CHROMOSOMES, maf=wildcards.maf, site=wildcards.site)
    else:
        return expand(rules.subset_angsd_maf.output, chrom=CHROMOSOMES, maf=wildcards.maf, site=wildcards.site)

def get_files_for_saf_estimation_byHabitat(wildcards):
    ref = rules.unzip_reference.output
    bams = expand(rules.create_bam_list_byHabitat.output, habitat = wildcards.habitat)
    return { 'bams' : bams, 'ref' : ref }

def get_habitat_saf_files(wildcards):
    all_saf_files = expand(rules.angsd_saf_likelihood_byHabitat.output.saf_idx, chrom=wildcards.chrom, habitat=HABITATS)
    first_hab = wildcards.hab_comb.split('_')[0]
    second_hab = wildcards.hab_comb.split('_')[1]
    saf1 = [x for x in all_saf_files if '_{0}_'.format(first_hab) in os.path.basename(x)]
    saf2 = [x for x in all_saf_files if '_{0}_'.format(second_hab) in os.path.basename(x)]
    return saf1 + saf2

def get_whatshap_phase_input(wildcards):
    ref = rules.unzip_reference.output
    vcf = rules.bcftools_split_samples.output
    bam = expand(rules.samtools_markdup.output.bam, sample = wildcards.sample)
    return { 'ref' : ref, 'vcf' : vcf, 'bam' : bam }
