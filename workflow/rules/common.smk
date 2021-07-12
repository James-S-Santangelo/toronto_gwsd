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

def get_representative_bam(wildcards):
    bam_index_files = expand(rules.index_bam.output, sample=SAMPLES)
    for i in bam_index_files:
        if REPRESENTATIVE_SAMPLE in i:
            bam = os.path.splitext(i)[0]
    return bam

def get_node_vcfs(wildcards):
    all_vcfs = expand(rules.bgzip_vcf.output, chrom=CHROMOSOMES, node=NODES)
    node_vcfs = [vcf for vcf in all_vcfs if wildcards.chrom in vcf]
    return node_vcfs

def get_node_tabix_files(wildcards):
    all_indices = expand(rules.tabix_node_vcf.output, chrom=CHROMOSOMES, node=NODES)
    node_indices = [i for i in all_indices if wildcards.chrom in i]
    return node_indices

def get_bed_to_subset(wildcards):
    all_bed_files = rules.get_fourfold_zerofold.output
    bed = [bed for bed in all_bed_files if wildcards.site in os.path.basename(bed)]
    return bed

def get_bams_for_angsd_gls(wildcards):
    if wildcards.sample_set == 'highQualSamples':
        return rules.create_bam_list_highQualSamples.output
    elif wildcards.sample_set == 'finalSamples_relatedRemoved':
        return rules.create_bam_list_allFinalSamples.output

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
    if wildcards.site == '0fold':
        return expand(rules.angsd_diversity_neutrality_stats.output, chrom=CHROMOSOMES, site='0fold', sample_set='finalSamples_relatedRemoved')
    elif wildcards.site == '4fold':
        return expand(rules.angsd_diversity_neutrality_stats.output, chrom=CHROMOSOMES, site='4fold', sample_set='finalSamples_relatedRemoved')
    else:
        return expand(rules.angsd_diversity_neutrality_stats.output, chrom=CHROMOSOMES, site='allSites', sample_set='finalSamples_relatedRemoved')

def get_angsd_sfs_toConcat(wildcards):
    if wildcards.site == '0fold':
        return expand(rules.angsd_estimate_sfs.output, chrom=CHROMOSOMES, site='0fold', sample_set='finalSamples_relatedRemoved')
    elif wildcards.site == '4fold':
        return expand(rules.angsd_estimate_sfs.output, chrom=CHROMOSOMES, site='4fold', sample_set='finalSamples_relatedRemoved')
    else:
        return expand(rules.angsd_estimate_sfs.output, chrom=CHROMOSOMES, site='allSites', sample_set='finalSamples_relatedRemoved')

def get_sites_for_angsd_index(wildcards):
    if wildcards.site == 'allSites':
        return rules.extract_angsd_allSites.output
    elif wildcards.site == '0fold' or wildcards.site == '4fold':
        return rules.split_angsd_sites_byChrom.output

def get_angsd_gl_toConcat(wildcards):
    if wildcards.site == '0fold' and wildcards.maf == '0.05' and wildcards.sample_set == 'highQualSamples':
        return expand(rules.subset_angsd_gl.output, site='0fold', maf='0.05', chrom=CHROMOSOMES, sample_set='highQualSamples')
    elif wildcards.site == '0fold' and wildcards.maf == '0.05' and wildcards.sample_set == 'finalSamples_relatedRemoved':
        return expand(rules.subset_angsd_gl.output, site='0fold', maf='0.05', chrom=CHROMOSOMES, sample_set='finalSamples_relatedRemoved')
    elif wildcards.site == '4fold' and wildcards.maf == '0.05' and wildcards.sample_set == 'highQualSamples':
        return expand(rules.subset_angsd_gl.output, site='4fold', maf='0.05', chrom=CHROMOSOMES, sample_set='highQualSamples')
    elif wildcards.site == '4fold' and wildcards.maf == '0.05' and wildcards.sample_set == 'finalSamples_relatedRemoved':
        return expand(rules.subset_angsd_gl.output, site='4fold', maf='0.05', chrom=CHROMOSOMES, sample_set='finalSamples_relatedRemoved')
    elif wildcards.site == 'allSites' and wildcards.maf == '0.05' and wildcards.sample_set == 'highQualSamples':
        return expand(rules.angsd_gl_allSites.output.gls, maf='0.05', chrom=CHROMOSOMES, sample_set='highQualSamples')
    elif wildcards.site == 'allSites' and wildcards.maf == '0.05' and wildcards.sample_set == 'finalSamples_relatedRemoved':
        return expand(rules.angsd_gl_allSites.output.gls, maf='0.05', chrom=CHROMOSOMES, sample_set='finalSamples_relatedRemoved')

def get_angsd_maf_toConcat(wildcards):
    if wildcards.site == '0fold' and wildcards.maf == '0.05' and wildcards.sample_set == 'highQualSamples':
        return expand(rules.subset_angsd_maf.output, site='0fold', maf='0.05', chrom=CHROMOSOMES, sample_set='highQualSamples')
    elif wildcards.site == '0fold' and wildcards.maf == '0.05' and wildcards.sample_set == 'finalSamples_relatedRemoved':
        return expand(rules.subset_angsd_maf.output, site='0fold', maf='0.05', chrom=CHROMOSOMES, sample_set='finalSamples_relatedRemoved')
    elif wildcards.site == '4fold' and wildcards.maf == '0.05' and wildcards.sample_set == 'highQualSamples':
        return expand(rules.subset_angsd_maf.output, site='4fold', maf='0.05', chrom=CHROMOSOMES, sample_set='highQualSamples')
    elif wildcards.site == '4fold' and wildcards.maf == '0.05' and wildcards.sample_set == 'finalSamples_relatedRemoved':
        return expand(rules.subset_angsd_maf.output, site='4fold', maf='0.05', chrom=CHROMOSOMES, sample_set='finalSamples_relatedRemoved')
    elif wildcards.site == 'allSites' and wildcards.maf == '0.05' and wildcards.sample_set == 'highQualSamples':
        return expand(rules.angsd_gl_allSites.output.mafs, maf='0.05', chrom=CHROMOSOMES, sample_set='highQualSamples')
    elif wildcards.site == 'allSites' and wildcards.maf == '0.05' and wildcards.sample_set == 'finalSamples_relatedRemoved':
        return expand(rules.angsd_gl_allSites.output.mafs, maf='0.05', chrom=CHROMOSOMES, sample_set='finalSamples_relatedRemoved')
