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

def get_fastas_to_concat(wildcards):
    if wildcards.gene == 'rbcl':
        return expand(rules.chloroplast_gene_consensus.output, sample=SAMPLES, gene='rbcl')
    elif wildcards.gene == 'matk':
        return expand(rules.chloroplast_gene_consensus.output, sample=SAMPLES, gene='matk')

def aggregate_ncbi_input(wildcards):
    checkpoint_output = checkpoints.download_nt_database.get(**wildcards).output[0]
    DB = glob_wildcards(os.path.join(checkpoint_output, 'nt.{db}.tar.gz')).db
    return expand('{0}/ncbi_nt_database/nt.{{db}}.{{ext}}'.format(PROGRAM_RESOURCE_DIR), db=DB, ext=['nhd', 'nhi', 'nhr', 'nin', 'nnd', 'nni', 'nog', 'nsq'])

def get_subset_bams_degeneracy_input(wildcards):
    """
    Returns the correct GLUE or Toronto BAM file
    """
    all_degen_bed_files = expand(rules.get_fourfold_zerofold.output, site=['0fold', '4fold'])
    regions = [x for x in all_degen_bed_files if wildcards.site in os.path.basename(x)]
    bam = expand(rules.samtools_markdup.output.bam, sample=wildcards.sample)
    idx = expand(rules.index_bam.output, sample = wildcards.sample)
    return { 'bam' : bam, 'idx' : idx, 'regions' : regions }

def get_vcfs_by_chrom(wildcards):
    return expand(rules.freebayes_call_variants.output, chrom=wildcards.chrom, i=FREEBAYES_CHUNKS)

def get_bed_to_subset(wildcards):
    all_bed_files = rules.get_fourfold_zerofold.output
    bed = [bed for bed in all_bed_files if wildcards.site in os.path.basename(bed)]
    return bed

def get_files_for_saf_estimation_byHabitat(wildcards):
    ref = rules.unzip_reference.output
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

def get_whatshap_phase_input(wildcards):
    ref = rules.unzip_reference.output
    vcf = rules.bcftools_split_samples.output
    bam = expand(rules.samtools_markdup.output.bam, sample = wildcards.sample)
    return { 'ref' : ref, 'vcf' : vcf, 'bam' : bam }

def get_dadi_sfs_input_files(wildcards):
    hab1 = wildcards.hab_comb.split('_')[0]
    hab2 = wildcards.hab_comb.split('_')[1]
    saf_files = expand(rules.angsd_saf_likelihood_byHabitat.output.saf_idx, habitat=HABITATS, site='4fold')  
    sfs_files = expand(rules.angsd_estimate_sfs_byHabitat.output, habitat=HABITATS, site='4fold') 
    saf_urban = [x for x in saf_files if '{0}'.format(hab1) in os.path.basename(x)]
    saf_rural = [x for x in saf_files if '{0}'.format(hab2) in os.path.basename(x)]
    sfs_urban = [x for x in sfs_files if '{0}'.format(hab1) in os.path.basename(x)]
    sfs_rural = [x for x in sfs_files if '{0}'.format(hab2) in os.path.basename(x)]
    ref = rules.unzip_reference.output
    return { 'saf_urban' : saf_urban , 'saf_rural' : saf_rural, 'sfs_urban' : sfs_urban, 'sfs_rural' : sfs_rural, 'ref' : ref }

def selscan_xpnsl_input(wildcards):
    if 'Urban' in wildcards.hab_comb:
        vcf = expand(rules.bcftools_splitVCF_byHabitat.output.vcf, chrom=wildcards.chrom, habitat='Urban')
    elif 'Suburban' in wildcards.hab_comb:
        vcf = expand(rules.bcftools_splitVCF_byHabitat.output.vcf, chrom=wildcards.chrom, habitat='Suburban')
    vcf_ref = expand(rules.bcftools_splitVCF_byHabitat.output.vcf, chrom=wildcards.chrom, habitat='Rural')
    return { 'vcf' : vcf, 'vcf_ref' : vcf_ref } 

def xpclr_input(wildcards):
    if 'Urban' in wildcards.hab_comb:
        pop1s = expand(rules.samples_byHabitat.output, chrom=wildcards.chrom, habitat='Urban')
    elif 'Suburban' in wildcards.hab_comb:
        pop1s = expand(rules.samples_byHabitat.output, chrom=wildcards.chrom, habitat='Suburban')
    pop2s = expand(rules.samples_byHabitat.output, chrom=wildcards.chrom, habitat='Rural')
    vcf = expand(rules.shapeit_phase.output.vcf, chrom=wildcards.chrom)
    return { 'vcf' : vcf, 'pop1s' : pop1s, 'pop2s' : pop2s } 


