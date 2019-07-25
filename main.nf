#!/usr/bin/env nextflow

/*
 * Authors:
 *      Mamana Mbiyavanga
 *
 *  On behalf of the H3ABionet Consortium
 *  2017
 *
 *
 * Description  : Create a test data from a vcf file
 * How to run: nextflow run main.nf --source file.vcf
 *
*/

// Show help emssage
//if (params.help){
//    helpMessage()
//    exit 0
//}

// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
nf_required_version = '19.04.0'
try {
    if( ! nextflow.version.matches(">= $nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
            "  Nextflow version $nf_required_version required! You are running v$workflow.nextflow.version.\n" +
            "  Pipeline execution will continue, but things may break.\n" +
            "  Please run `nextflow self-update` to update Nextflow.\n" +
            "============================================================"
}

//source  = Channel.from(params.source)
source  = Channel.fromPath(params.source_dataset).view()
subset_size = params.subset_size

process get_snp_list {
    tag "get_snp_list_${base}_${subset_size}"
//    publishDir "${params.outDir}/testdata", overwrite: true, mode:'copy', pattern: "*map"
    input:
        file(vcf_file) from source
    output:
        set file(vcf_file), file(subset_map_ref), file(subset_map_target) into get_snp_list,get_snp_list_1,get_snp_list_2,get_snp_list_3
        file(chromFile) into chromosome
    script:
        base = file(vcf_file.baseName).baseName
        vcf_ref_chunk = "${base}_${params.chunk}.vcf.gz"
        chromFile = "${base}.chrom"
        subset_map_ref = "refPanel_testdata_${subset_size}.map"
        subset_bed_ref = "refPanel_testdata_${subset_size}.bed"
        subset_map_target = "target_testdata_${subset_size}"
        """
        tabix ${vcf_file}
        bcftools view --regions ${params.chunk} ${vcf_file} -Oz -o ${vcf_ref_chunk}
        bcftools query -f '%CHROM\\t%POS\\t%POS\\n' ${vcf_ref_chunk} > ${base}.map
        awk -F' ' '{print \$1}' ${base}.map | sort -n | uniq > ${chromFile}
        awk -F'\\t' '{print \$1"\\t"\$2-1"\\t"\$3"\\t"\$1" dna:chromosome chromosome:GRCh37:"\$1":"\$2-1":"\$3":1"}' ${base}.map > ${subset_map_ref}
        sort -R ${base}.map | tail -n ${subset_size} > ${subset_map_target}
        """
}

chromosomes = file(chromosome.toSortedList().val[0]).readLines().unique().collect { it as int }.sort()

process subset_vcf {
    tag "subset_vcf_${base}_${subset_size}"
    publishDir "${params.publishDir}", overwrite: true, mode:'copy', pattern: "${vcf_target}"
    input:
        set file(vcf_file), file(subset_map_ref), file(subset_map_target) from get_snp_list
    output:
        set file(vcf_file), file(subset_map_ref), file(subset_map_target), file(vcf_ref), file(vcf_target) into subset_vcf_ref,subset_vcf_ref_1
    script:
        base = file(vcf_file.baseName).baseName
        vcf_ref = "refPanel_testdata.vcf.gz"
        vcf_target = "target_testdata.vcf.gz"
        """
        tabix ${vcf_file}
        bcftools view --regions-file ${subset_map_ref} ${vcf_file} -Oz -o ${vcf_ref}
        tabix ${vcf_ref}
        bcftools view --regions-file ${subset_map_target} ${vcf_ref} -Oz -o ${vcf_target}
        """
}


process split_vcf_to_chrm {
    tag "split_${base}_${chrm}"
    label "medium"
    input:
        each chrm from chromosomes
        set file(vcf_file), file(subset_map_ref), file(subset_map_target), file(vcf_ref), file(vcf_target) from subset_vcf_ref
    output:
        set chrm, file(vcf_ref_chrm) into split_vcf_to_chrm
    script:
        base = file(vcf_ref.baseName).baseName
        vcf_ref_chrm = "${base}_${chrm}.vcf.gz"
        """
        tabix ${vcf_ref}
        bcftools view \
            --regions ${chrm} \
            -m2 -M2 -v snps \
            ${vcf_ref} \
            -Oz -o ${vcf_ref_chrm}
        """
}


process phase_vcf_chrm {
    tag "phase_${base}_${chrm}"
    publishDir "${params.publishDir}", overwrite: true, mode:'copy', pattern: "*.vcf.gz"

    input:
        set chrm, file(vcf_file_chrm) from split_vcf_to_chrm
    output:
        set chrm, file(vcf_file_chrm), file("${file_out}.vcf.gz") into phase_vcf_chrm
    script:
        base = file(vcf_file_chrm.baseName).baseName
        file_out = "${base}_phased"
        """
        nblines=\$(zcat ${vcf_file_chrm} | grep -v '^#' | wc -l)
        if (( \$nblines > 0 ))
        then
            tabix ${vcf_file_chrm}
            eagle \
                --vcfTarget=${vcf_file_chrm} \
                --geneticMapFile=${params.eagle_genetic_map} \
                --vcfRef=${params.reference_vcf} \
                --vcfOutFormat=z \
                --chrom=${chrm} \
                --outPrefix=${file_out} 2>&1 | tee ${file_out}.log
            if [ ! -f "${file_out}.vcf.gz" ]; then
                touch ${file_out}.vcf && bgzip -f ${file_out}.vcf
            fi
        else
            touch ${file_out}.vcf && bgzip -f ${file_out}.vcf
        fi
        """
}


process vcf_to_m3vcf {
    tag "m3vcf_${ref}_${chrm}"
    label "medmem"
    publishDir "${params.publishDir}", overwrite: true, mode:'copy', pattern: "${base}.m3vcf.gz"
    input:
        set chrm, file(vcf_chrm), file(vcf_chrm_phased) from phase_vcf_chrm
    output:
        set chrm, file(vcf_chrm), file(vcf_chrm_phased), file(m3vcf_chrm) into vcf_to_m3vcf
    script:
        base = file(vcf_chrm_phased.baseName).baseName
        m3vcf_chrm = "${base}.m3vcf.gz"
        """
        minimac3 \
            --refHaps ${vcf_chrm_phased} \
            --processReference \
            --prefix ${base}
        """
}


//"""
//Generate max and min for each chromosome from a bed file
//"""
//process get_coordinates {
//    tag "get_coordinates_${base}"
//    label "medmem"
//    publishDir "${params.outDir}/testdata", overwrite: true, mode:'copy'
//    input:
//        set file(vcf_file), file(subset_map_ref), file(subset_map_target) from get_snp_list_1
//    output:
//        set file(bed_output), file(subset_map_ref) into get_coordinates
//    script:
//        base = subset_map_target.baseName
//        bedFile = subset_map_target
//        bed_output = "${base}_coordinates.bed"
//        template "get_coordinates_from_bed.py"
//}


"""
Subset a fasta file based on a bed file
"""
process subset_fasta {
    tag "subset_fasta_${subset_map_ref.baseName}"
    label "medmem"

    input:
    set file(vcf_file), file(subset_map_ref), file(subset_map_target), file(vcf_ref), file(vcf_target) from subset_vcf_ref_1

    output:
    file(ref_fasta) into subset_fasta

    script:
    ref_fasta = "hg19_testdata_temp.fasta"
    """
    bedtools getfasta -fi ${params.reference_genome} -bed ${vcf_target} -fo ${ref_fasta}
    """
}


"""
Fix header of fasta file
"""
process fix_fasta {
    tag "fix_fasta_${fasta_in.baseName}"
    label "medmem"

    input:
    file(fasta_in) from subset_fasta

    output:
    file(fasta_out) into fix_fasta

    script:
    fasta_out = "hg19_testdata.fasta"
    template "fix_fasta_header.py"
}

"""
Index fasta file
"""
process index_fasta {
    tag "index_fasta_${fasta.baseName}"
    label "medmem"
    publishDir "${params.publishDir}", overwrite: true, mode:'copy', pattern: "${fasta}*"

    input:
    file(fasta) from fix_fasta

    output:
    set file(fasta), file("${fasta}.fai") into index_fasta

    script:
    """
    samtools faidx ${fasta}
    """
}


"""
Generate test genetic map 
"""
process subset_genetic_map {
    tag "subset_genetic_map_${base}"
    label "medmem"
    publishDir "${params.publishDir}", overwrite: true, mode:'copy', pattern: "${out_map}*"

    input:
    set file(vcf_file), file(map_ref), file(map_target) from get_snp_list_2

    output:
    set file(map_ref), file(out_map) into subset_genetic_map

    script:
    map = file(params.eagle_genetic_map)
    subset_map = map_ref
    base = file(file(map).baseName).baseName
    out_map = "${base}_testdata.txt.gz"
    template "subset_genetic_map.py"
}

"""
Generate test chip 
"""
process subset_chip {
    tag "subset_chip_${base}"
    label "medmem"
    publishDir "${params.publishDir}", overwrite: true, mode:'copy', pattern: "${out_chip}*"

    input:
    set file(vcf_file), file(map_ref), file(map_target) from get_snp_list_3

    output:
    set file(map_ref), file(out_chip) into subset_chip

    script:
    chip = file(params.h3achip)
    base = file(file(chip).baseName).baseName
    subset_map = map_ref
    out_chip = "${base}_testdata.csv"
    template "subset_chip.py"
}
