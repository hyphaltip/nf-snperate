#!/usr/bin/env nextflow

/*
========================================================================================
    SNP Calling Pipeline using GATK and bwa-mem2
    Translated from bash pipeline scripts
========================================================================================
    Github : https://github.com/jstajich/snperate
    Contact: Jason Stajich
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

log.info "========================================="
log.info "  SNPERATE - SNP Calling Pipeline"
log.info "========================================="
log.info "samples          : ${params.samples}"
log.info "population_yaml : ${params.population_yaml}"
log.info "outdir          : ${params.outdir}"
log.info "filter_method   : ${params.filter_method}"
log.info "genome_version : ${params.genome_version}"
log.info "========================================="

/*
========================================================================================
    PARSE SAMPLES CSV
========================================================================================
*/

def samples = []
file(params.samples).eachLine { line ->
    line = line.trim()
    if (line && !line.startsWith('#')) {
        def parts = line.split(',')
        if (parts.size() >= 2) {
            def sample_id = parts[0].trim()
            def fastq_paths = parts[1].trim().split(';')*.trim()
            def platform = parts.size() > 2 ? parts[2].trim() : "Illumina"
            samples << [id: sample_id, fastqs: fastq_paths, platform: platform]
        }
    }
}

log.info "Found ${samples.size()} samples"

/*
========================================================================================
    PARSE POPULATIONS YAML
========================================================================================
*/

def populations = [:]
def in_populations = false
def current_pop = null

file(params.population_yaml).eachLine { line ->
    line = line.trim()
    if (line.startsWith('Populations:')) {
        in_populations = true
    } else if (in_populations && line.endsWith(':')) {
        current_pop = line.replace(':', '').trim()
        populations[current_pop] = []
    } else if (in_populations && current_pop && line.startsWith('-')) {
        def strain = line.replace('-', '').trim()
        if (strain) {
            populations[current_pop] << strain
        }
    }
}

// Add "All" population
populations['All'] = samples.collect { it.id }

log.info "Found ${populations.size()} populations: ${populations.keySet().join(', ')}"

/*
========================================================================================
    CHANNELS
========================================================================================
*/

sample_ch = Channel.fromList(samples)

/*
========================================================================================
    PROCESS: INDEX REFERENCE GENOME
========================================================================================
*/

process INDEX_REFERENCE {
    label 'index'
    cpus 4
    memory '8 GB'
    publishDir "${params.genome_folder}", mode: 'copy', overwrite: false

    input:
    path refgenome

    output:
    path "*.bwt", emit: bwa_index
    path "*.fai", emit: fai
    path "*.dict", emit: dict
    path "chrom_nums.txt", emit: chroms

    script:
    def base = refgenome.simpleName
    """
    bwa index $refgenome
    samtools faidx $refgenome
    samtools dict $refgenome > ${base}.dict
    grep '>' $refgenome | sed 's/>//' | cut -d' ' -f1 > chrom_nums.txt
    """
}

/*
========================================================================================
    PROCESS: ALIGN READS WITH BWA-MEM2
========================================================================================
*/

process ALIGN_READS {
    label 'align'
    cpus 8
    memory '32 GB'
    stageInMode 'copy'
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(sample_id), path(fastqs)
    path refgenome

    output:
    tuple val(sample_id), path("${sample_id}.cram"), emit: crams
    path "${sample_id}.cram.crai", emit: crais

    script:
    def rg_line = "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tLB:${sample_id}\\tPL:${params.rg_platform}\\tCN:${params.rg_center}"
    def fastq_input = fastqs.size() >= 2 ? "${fastqs[0]} ${fastqs[1]}" : "${fastqs[0]}"

    """
    bwa-mem2 mem -t ${task.cpus} -R "${rg_line}" $refgenome $fastq_input | \\
        samtools sort -@ ${task.cpus} -O cram -T tmp_${sample_id} -o ${sample_id}.cram
    samtools index -@ ${task.cpus} ${sample_id}.cram
    """
}

/*
========================================================================================
    PROCESS: CALL GVCF WITH GATK HAPLOTYPECALLER
========================================================================================
*/

process CALL_GVCF {
    label 'gvcf'
    cpus 8
    memory '48 GB'
    stageInMode 'copy'

    input:
    tuple val(sample_id), path(cram), path(crai)
    path refgenome

    output:
    tuple val(sample_id), path("${sample_id}.g.vcf.gz"), emit: gvcfs
    path "${sample_id}.g.vcf.gz.tbi", emit: tbis

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" HaplotypeCaller \\
        --emit-ref-confidence GVCF \\
        --sample-ploidy 1 \\
        --input ${sample_id}.cram \\
        --reference ${refgenome} \\
        --output ${sample_id}.g.vcf.gz \\
        --native-pair-hmm-threads ${task.cpus} \\
        --sample-name ${sample_id} \\
        -G StandardAnnotation \\
        -G AS_StandardAnnotation \\
        -G StandardHCAnnotation
    """
}

/*
========================================================================================
    PROCESS: IMPORT GVCFS TO GENOMICSDB FOR JOINT GENOTYPING
========================================================================================
*/

process IMPORT_GVCF {
    label 'genomicsdb'
    cpus 8
    memory '32 GB'
    stageInMode 'copy'

    input:
    path gvcfs
    path tbis
    val interval
    val pop_name

    output:
    path "genomicsdb", emit: genomicsdb

    script:
    def gvcfs_list = gvcfs.collect { "-V $it" }.join(' ')
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" GenomicsDBImport \\
        --genomicsdb-workspace-path genomicsdb \\
        $gvcfs_list \\
        --intervals ${interval} \\
        --reader-threads ${task.cpus} \\
        --tmp-dir ./
    """
}

/*
========================================================================================
    PROCESS: JOINT GENOTYPING WITH GATK
========================================================================================
*/

process JOINT_GENOTYPE {
    label 'genotype'
    cpus 4
    memory '24 GB'

    input:
    path genomicsdb
    path refgenome
    val interval
    val pop_name

    output:
    path "joint_${interval}_${pop_name}.vcf.gz", emit: joint_vcf

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" GenotypeGVCFs \\
        --reference ${refgenome} \\
        --output joint_${interval}_${pop_name}.vcf.gz \\
        -V gendb://${genomicsdb} \\
        --tmp-dir ./
    """
}

/*
========================================================================================
    PROCESS: SELECT SNP VARIANTS
========================================================================================
*/

process SELECT_SNP {
    label 'select'
    cpus 2
    memory '8 GB'

    input:
    path vcf
    path refgenome
    val interval
    val pop_name

    output:
    path "snp_${interval}_${pop_name}.vcf.gz", emit: snp_vcf

    script:
    """
    gatk SelectVariants \\
        -R ${refgenome} \\
        --variant ${vcf} \\
        -O snp_${interval}_${pop_name}.vcf.gz \\
        --restrict-alleles-to BIALLELIC \\
        --select-type-to-include SNP \\
        --create-output-variant-index false
    """
}

/*
========================================================================================
    PROCESS: FILTER SNP VARIANTS WITH GATK
========================================================================================
*/

process FILTER_SNP_GATK {
    label 'filter'
    cpus 2
    memory '8 GB'

    input:
    path vcf
    path refgenome
    val interval
    val pop_name

    output:
    path "filter_snp_${interval}_${pop_name}.vcf.gz", emit: filtered_snp

    script:
    """
    gatk VariantFiltration \\
        --variant ${vcf} \\
        -R ${refgenome} \\
        --output filter_snp_${interval}_${pop_name}.vcf.gz \\
        --cluster-window-size 10 \\
        --filter-expression "QD < 2.0" --filter-name QualByDepth \\
        --filter-expression "MQ < 40.0" --filter-name MapQual \\
        --filter-expression "QUAL < 100" --filter-name QScore \\
        --filter-expression "SOR > 4.0" --filter-name StrandOddsRatio \\
        --filter-expression "FS > 60.0" --filter-name FisherStrandBias \\
        --missing-values-evaluate-as-failing \\
        --create-output-variant-index false
    """
}

/*
========================================================================================
    PROCESS: FILTER SNP VARIANTS WITH BCFTOOLS
========================================================================================
*/

process FILTER_SNP_BCFTOOLS {
    label 'filter'
    cpus 4
    memory '8 GB'

    input:
    path vcf
    val interval
    val pop_name

    output:
    path "filter_bcf_snp_${interval}_${pop_name}.vcf.gz", emit: filtered_snp

    script:
    """
    bcftools view -i '(QD < 2.0) || (MQ < 40.0) || (QUAL < 100) || (SOR > 4.0) || (FS > 60.0)' \\
        -Oz -o filter_bcf_snp_${interval}_${pop_name}.vcf.gz ${vcf}
    tabix filter_bcf_snp_${interval}_${pop_name}.vcf.gz
    """
}

/*
========================================================================================
    PROCESS: SELECT PASSING VARIANTS
========================================================================================
*/

process SELECT_PASS {
    label 'filter'
    cpus 2
    memory '4 GB'

    input:
    path vcf
    path refgenome
    val interval
    val pop_name

    output:
    path "pass_${interval}_${pop_name}.vcf.gz", emit: pass_vcf

    script:
    """
    gatk SelectVariants \\
        -R ${refgenome} \\
        --variant ${vcf} \\
        --output pass_${interval}_${pop_name}.vcf.gz \\
        --exclude-filtered \\
        --create-output-variant-index false
    """
}

/*
========================================================================================
    PROCESS: COMBINE SLICED VCFs
========================================================================================
*/

process COMBINE_VCF {
    label 'vcf'
    cpus 4
    memory '8 GB'
    publishDir "${params.outdir}/vcf", mode: 'copy'

    input:
    path vcfs
    val pop_name

    output:
    path "combined_${pop_name}.vcf.gz", emit: combined_vcf
    path "combined_${pop_name}.vcf.gz.tbi", emit: combined_tbi

    script:
    """
    bcftools concat -a -Oz -o combined_${pop_name}.vcf.gz $vcfs
    tabix combined_${pop_name}.vcf.gz
    """
}

/*
========================================================================================
    PROCESS: CNV DETECTION WITH MOSDEPTH
========================================================================================
*/

process MOSDEPTH {
    label 'cnv'
    cpus 4
    memory '16 GB'
    publishDir "${params.outdir}/coverage", mode: 'copy'

    input:
    tuple val(sample_id), path(cram), path(crai)
    path refgenome
    val window

    output:
    path "mosdepth/${sample_id}.${window}bp.*", emit: coverage

    script:
    """
    mosdepth -f ${refgenome} -T 1,10,50,100,200 -n --by ${window} -t ${task.cpus} ${sample_id}.${window}bp ${cram}
    """
}

/*
========================================================================================
    WORKFLOW DEFINITION
========================================================================================
*/

workflow SNPERATE_PIPELINE {

    // Reference genome file
    def genome_file = "${params.genome_folder}/${params.genome_name}_${params.genome_version}_Genome.fasta"
    def ref_file = file(genome_file)

    if (!ref_file.exists()) {
        log.error "Reference genome not found at: ${genome_file}"
        log.error "Please set correct genome path in nextflow.config"
        exit 1
    }

    // Create ref genome channel
    ref_ch = Channel.value(ref_file)

    // Check for existing index files
    def genome_dir = file(params.genome_folder)
    def has_index = genome_dir.exists() && genome_dir.listFiles().findAll { it.name.endsWith('.bwt') }.size() > 0

    if (!has_index) {
        INDEX_REFERENCE(ref_ch)
        ref_fai = INDEX_REFERENCE.out.fai
        ref_dict = INDEX_REFERENCE.out.dict
        chrom_ch = INDEX_REFERENCE.out.chroms
    } else {
        ref_fai = Channel.fromPath("${params.genome_folder}/*.fai").first()
        ref_dict = Channel.fromPath("${params.genome_folder}/*.dict").first()
        chrom_ch = Channel.fromPath("${params.genome_folder}/chrom_nums.txt").first()
    }

    // Align reads
    ALIGN_READS(sample_ch, ref_ch)

    // Call GVCFs
    cram_crai_ch = ALIGN_READS.out.crams.join(ALIGN_READS.out.crais)
    CALL_GVCF(cram_crai_ch, ref_ch)

    // Get chromosome list
    chrom_list = chrom_ch.splitText().map { it.trim() }.findAll { it }

    // Create interval channel from chromosomes
    interval_ch = Channel.from(chrom_list)

    // Create population channel
    pop_ch = Channel.from(populations.collect { [it.key, it.value] })

    // Group GVCFs by population
    gvcf_grouped = CALL_GVCF.out.gvcfs
        .combine(pop_ch) { sample_gvcf, pop ->
            [sample_gvcf[0], sample_gvcf[1], pop[0], pop[1]]
        }
        .filter { sample_id, gvcf, pop_name, strains -> sample_id in strains }
        .map { sample_id, gvcf, pop_name, strains -> [pop_name, gvcf] }
        .groupTuple()

    tbi_grouped = CALL_GVCF.out.tbis
        .combine(pop_ch) { sample_tbi, pop ->
            [sample_tbi[0], sample_tbi[1], pop[0], pop[1]]
        }
        .filter { sample_id, tbi, pop_name, strains -> sample_id in strains }
        .map { sample_id, tbi, pop_name, strains -> [pop_name, tbi] }
        .groupTuple()

    // Create input for import: [pop_name, [gvcfs], [tbis], interval]
    import_input = gvcf_grouped
        .combine(tbi_grouped) { pop1, gvcfs, pop2, tbis ->
            [pop1, gvcfs, tbis]
        }
        .filter { pop_name, gvcfs, tbis, pop_name2, gvcfs2, tbis2 -> pop_name == pop_name2 }
        .map { pop_name, gvcfs, tbis, pop_name2, gvcfs2, tbis2 ->
            [pop_name, gvcfs, tbis]
        }
        .combine(interval_ch) { pop_name, gvcfs, tbis, interval ->
            [pop_name, gvcfs, tbis, interval]
        }

    // Import to GenomicsDB
    IMPORT_GVCF(import_input.map { it[1] }, import_input.map { it[2] }, import_input.map { it[3] }, import_input.map { it[0] })

    // Joint genotyping
    joint_input = IMPORT_GVCF.out
        .combine(ref_ch)
        .combine(interval_ch)
        .combine(pop_ch.map { it[0] })

    JOINT_GENOTYPE(joint_input.map { it[0] }, joint_input.map { it[1] }, joint_input.map { it[2] }, joint_input.map { it[3] })

    // Select SNPs
    select_input = JOINT_GENOTYPE.out
        .combine(ref_ch)
        .combine(interval_ch)
        .combine(pop_ch.map { it[0] })

    SELECT_SNP(select_input.map { it[0] }, select_input.map { it[1] }, select_input.map { it[2] }, select_input.map { it[3] })

    // Filter SNPs
    filter_input = SELECT_SNP.out
        .combine(ref_ch)
        .combine(interval_ch)
        .combine(pop_ch.map { it[0] })

    if (params.filter_method == "gatk") {
        FILTER_SNP_GATK(filter_input.map { it[0] }, filter_input.map { it[1] }, filter_input.map { it[2] }, filter_input.map { it[3] })
        filtered_ch = FILTER_SNP_GATK.out.filtered_snp
    } else {
        FILTER_SNP_BCFTOOLS(filter_input.map { it[0] }, filter_input.map { it[2] }, filter_input.map { it[3] })
        filtered_ch = FILTER_SNP_BCFTOOLS.out.filtered_snp
    }

    // Select passing variants
    pass_input = filtered_ch
        .combine(ref_ch)
        .combine(interval_ch)
        .combine(pop_ch.map { it[0] })

    SELECT_PASS(pass_input.map { it[0] }, pass_input.map { it[1] }, pass_input.map { it[2] }, pass_input.map { it[3] })

    // Combine VCFs per population
    pass_by_pop = SELECT_PASS.out.pass_vcf
        .map { vcf -> [vcf.name.split('_')[2], vcf] }
        .groupTuple()

    COMBINE_VCF(pass_by_pop.map { it[1] }, pass_by_pop.map { it[0] })

    // CNV detection
    window_ch = Channel.fromList(params.mosdepth_windows)
    MOSDEPTH(cram_crai_ch.combine(ref_ch).combine(window_ch))
}

/*
========================================================================================
    ENTRY POINT
========================================================================================
*/

workflow {
    SNPERATE_PIPELINE()
}
