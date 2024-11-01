# Python functions used throughout snakemake workflow

def run_fast_scandir(dir, sample):
    """
    Helper function to perform fast recursive search
    """
    subfolders, files = [], []

    for f in os.scandir(dir):
        if f.is_dir():
            subfolders.append(f.path)
        if f.is_file():
            if re.findall(r'^{0}_'.format(sample), f.name):
                files.append(f.path)
    for dir in list(subfolders):
        sf, f = run_fast_scandir(dir, sample)
        subfolders.extend(sf)
        files.extend(f) 
    return subfolders, sorted(files)

def get_raw_reads(wildcards):
    """
    Recursively search for forward and reverse reads for sample
    """
    folders, reads = run_fast_scandir(RAW_READ_DIR, wildcards.sample)
    R1 = reads[0]
    R2 = reads[1]
    return { 'read1' : R1, 'read2' : R2  }

def get_subset_bams_degeneracy_input(wildcards):
    """
    Returns the correct GLUE or Toronto BAM file
    """
    all_degen_bed_files = expand(rules.create_bed_from_degenotate.output, site=['0fold', '4fold'])
    regions = [x for x in all_degen_bed_files if wildcards.site in os.path.basename(x)]
    bam = expand(rules.samtools_markdup.output.bam, sample=wildcards.sample)
    idx = expand(rules.index_bam.output, sample = wildcards.sample)
    return { 'bam' : bam, 'idx' : idx, 'regions' : regions }

def get_vcfs_by_chrom(wildcards):
    return expand(rules.freebayes_call_variants.output, chrom=wildcards.chrom, i=FREEBAYES_CHUNKS)

def get_bed_to_subset(wildcards):
    bed = expand(rules.create_bed_from_degenotate.output, site=wildcards.site)
    return bed

def get_files_for_saf_estimation_byHabitat(wildcards):
    ref = REFERENCE_GENOME
    bams = expand(rules.create_bam_list_byHabitat.output, habitat = wildcards.habitat, site=wildcards.site)
    sites = rules.convert_sites_for_angsd.output
    idx = rules.angsd_index_degenerate_sites.output
    chroms = config['chromosomes']
    return { 'bams' : bams, 'ref' : ref, 'sites' : sites, 'idx' : idx, 'chroms' : chroms }

def get_habitat_saf_files(wildcards):
    all_saf_files = expand(rules.angsd_saf_likelihood_byHabitat.output.saf_idx, habitat=HABITATS, site=wildcards.site)
    first_hab = wildcards.hab_comb.split('_')[0]
    second_hab = wildcards.hab_comb.split('_')[1]
    saf1 = [x for x in all_saf_files if '{0}'.format(first_hab) in os.path.basename(x)]
    saf2 = [x for x in all_saf_files if '{0}'.format(second_hab) in os.path.basename(x)]
    return saf1 + saf2

def get_habitat_saf_files_allSites(wildcards):
    all_saf_files = expand(rules.angsd_saf_likelihood_byHabitat_allSites.output.saf_idx, habitat=HABITATS, chrom=wildcards.chrom)
    first_hab = wildcards.hab_comb.split('_')[0]
    second_hab = wildcards.hab_comb.split('_')[1]
    saf1 = [x for x in all_saf_files if '{0}'.format(first_hab) in os.path.basename(x)]
    saf2 = [x for x in all_saf_files if '{0}'.format(second_hab) in os.path.basename(x)]
    return saf1 + saf2

def get_files_for_saf_estimation_byPopulation(wildcards):
    ref = REFERENCE_GENOME
    bams = expand(rules.create_bam_list_byPop_multiInd.output, popu=wildcards.popu, site=wildcards.site)
    sites = rules.convert_sites_for_angsd.output
    idx = rules.angsd_index_degenerate_sites.output
    chroms = config['chromosomes']
    return { 'bams' : bams, 'ref' : ref, 'sites' : sites, 'idx' : idx, 'chroms' : chroms }

def get_population_saf_files(wildcards):
    all_saf_files = expand(rules.angsd_saf_likelihood_byPopulation.output.saf_idx, popu=POPS_MULTI_IND, site='4fold')
    pop1 = wildcards.pop_comb.split('_')[0]
    pop2 = wildcards.pop_comb.split('_')[1]
    saf1 = [x for x in all_saf_files if os.path.basename(x).startswith('{0}_'.format(pop1))]
    saf2 = [x for x in all_saf_files if os.path.basename(x).startswith('{0}_'.format(pop2))]
    return saf1 + saf2 

def get_whatshap_phase_input(wildcards):
    ref = REFERENCE_GENOME
    vcf = rules.bcftools_split_samples.output
    bam = expand(rules.samtools_markdup.output.bam, sample = wildcards.sample)
    return { 'ref' : ref, 'vcf' : vcf, 'bam' : bam }

def selscan_xpnsl_input(wildcards):
    if wildcards.hab_comb == 'Urban_Rural':
        vcf = expand(rules.bcftools_splitVCF_byHabitat.output.vcf, chrom=wildcards.chrom, habitat='Urban')
        vcf_ref = expand(rules.bcftools_splitVCF_byHabitat.output.vcf, chrom=wildcards.chrom, habitat='Rural')
    elif wildcards.hab_comb == 'Urban_Suburban':
        vcf = expand(rules.bcftools_splitVCF_byHabitat.output.vcf, chrom=wildcards.chrom, habitat='Urban')
        vcf_ref = expand(rules.bcftools_splitVCF_byHabitat.output.vcf, chrom=wildcards.chrom, habitat='Suburban')
    elif wildcards.hab_comb == 'Suburban_Rural':
        vcf = expand(rules.bcftools_splitVCF_byHabitat.output.vcf, chrom=wildcards.chrom, habitat='Suburban')
        vcf_ref = expand(rules.bcftools_splitVCF_byHabitat.output.vcf, chrom=wildcards.chrom, habitat='Rural')
    return { 'vcf' : vcf, 'vcf_ref' : vcf_ref } 

def bcftools_splitVCF_permuted_input(wildcards):
    if wildcards.habitat == "Urban":
       samples_allChroms = expand(rules.permuted_samples_byHabitat.output.urb, n=wildcards.n)
       samples = [x for x in samples_allChroms if wildcards.chrom in x]
    if wildcards.habitat == "Rural":
       samples_allChroms = expand(rules.permuted_samples_byHabitat.output.rur, n=wildcards.n)
       samples = [x for x in samples_allChroms if wildcards.chrom in x]
    return(samples)

def get_windowed_xpnsl_input_files(wildcards):
    if wildcards.hab_comb == 'Urban_Rural':
        norm = expand(rules.norm_xpnsl.output, chrom=CHROMOSOMES, hab_comb='Urban_Rural')
    elif wildcards.hab_comb == 'Urban_Suburban':
        norm = expand(rules.norm_xpnsl.output, chrom=CHROMOSOMES, hab_comb='Urban_Suburban')
    elif wildcards.hab_comb == 'Suburban_Rural':
        norm = expand(rules.norm_xpnsl.output, chrom=CHROMOSOMES, hab_comb='Suburban_Rural')
    return norm

def get_windowed_singPop_hapstats_input_files(wildcards):
    if wildcards.stat == "xpnsl":
        norm = expand(rules.norm_xpnsl.output, chrom=CHROMOSOMES, hab_comb='Urban_Rural')
    elif wildcards.stat == "nsl":
        norm = expand(rules.norm_nsl.output, chrom=CHROMOSOMES, habitat=['Urban', 'Rural'])
    elif wildcards.stat == "ihs":
        norm = expand(rules.norm_ihs.output, chrom=CHROMOSOMES, habitat=['Urban', 'Rural'])
    elif wildcards.stat == "ihh12":
        norm = expand(rules.norm_ihh_OneTwo.output, chrom=CHROMOSOMES, habitat=['Urban', 'Rural'])
    return norm

