/*
========================================================================================
    WGS analysis pipeline: fastp -> bwa-mem2 -> markdup -> deepvariant -> glnexus
========================================================================================
*/

// Use DSL2 for modern syntax
nextflow.enable.dsl=2

// -- Log pipeline parameters
log.info """
         W G S - A N A L Y S I S   P I P E L I N E (v2)
         =============================================
         Input Samplesheet : ${params.input}
         Reference FASTA   : ${params.ref_fasta}
         BWA Index Prefix  : ${params.bwa_index}
         Output Directory  : ${params.outdir}
         ---
         Run Joint Calling : ${params.joint_calling}
         Profile           : ${workflow.profile}
         =============================================
         """
         .stripIndent()

// -- Log container versions
log.info """
         C O N T A I N E R   V E R S I O N S
         =============================================
         FASTP             : ${params.fastp_container}
         BWA-MEM2          : ${params.bwa_mem2_container}
         SAMTOOLS          : ${params.samtools_container}
         DEEPVARIANT       : ${params.deepvariant_container}
         GLNEXUS           : ${params.glnexus_container}
         MULTIQC           : ${params.multiqc_container}
         =============================================
         """
         .stripIndent()

/*
 * ====================================================================================
 *   Define Processes
 * ====================================================================================
 */

process FASTP {
    tag "FASTP on ${sample_id}"
    publishDir "${params.outdir}/fastp/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.trimmed.R1.fastq.gz"), path("${sample_id}.trimmed.R2.fastq.gz"), emit: trimmed_reads
    path "${sample_id}.fastp.html", emit: html_report
    path "${sample_id}.fastp.json", emit: json_report

    script:
    """
    fastp \\
        -i ${reads[0]} \\
        -I ${reads[1]} \\
        -o ${sample_id}.trimmed.R1.fastq.gz \\
        -O ${sample_id}.trimmed.R2.fastq.gz \\
        -h ${sample_id}.fastp.html \\
        -j ${sample_id}.fastp.json \\
        -w ${task.cpus}
    """
}

process ALIGN_BWA_MEM2 {
    tag "BWA-MEM2 on ${sample_id}"
    publishDir "${params.outdir}/bam/${sample_id}", mode: 'copy', pattern: '*.{flagstat.txt,sorted.bam,sorted.bam.bai}'

    input:
    tuple val(sample_id), path(reads)
    path ref_fasta
    path bwa_index_prefix

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: bam_bai
    path "${sample_id}.flagstat.txt", emit: flagstat

    script:
    def bwa_index = bwa_index_prefix.toString()
    def rg_string = "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA"

    """
    bwa-mem2 mem \\
        -t ${task.cpus} \\
        -R '${rg_string}' \\
        ${bwa_index} \\
        ${reads[0]} \\
        ${reads[1]} | \\
    samtools sort -@ ${task.cpus-1} -o ${sample_id}.sorted.bam -

    samtools index ${sample_id}.sorted.bam
    samtools flagstat ${sample_id}.sorted.bam > ${sample_id}.flagstat.txt
    """
}

process MARK_DUPLICATES {
    tag "Mark Duplicates on ${sample_id}"
    publishDir "${params.outdir}/bam/${sample_id}", mode: 'copy', pattern: '*.dedup.bam*'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.dedup.bam"), path("${sample_id}.dedup.bam.bai"), emit: dedup_bam_bai

    script:
    """
    samtools markdup \\
        -@ ${task.cpus} \\
        -r \\
        ${bam} \\
        ${sample_id}.dedup.bam

    samtools index ${sample_id}.dedup.bam
    """
}

process DEEPVARIANT {
    tag "DeepVariant on ${sample_id}"
    publishDir "${params.outdir}/variants/single_sample/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path ref_fasta

    output:
    tuple val(sample_id), path("*.vcf.gz"), path("*.gvcf.gz"), path("*.vcf.gz.tbi"), path("*.html"), emit: results

    script:
    def ref_fasta_name = ref_fasta.getName()
    """
    ln -s ${ref_fasta} ${ref_fasta_name}
    samtools faidx ${ref_fasta_name}

    /opt/deepvariant/bin/run_deepvariant \\
        --model_type=WGS \\
        --ref=${ref_fasta_name} \\
        --reads=${bam} \\
        --output_vcf=${sample_id}.deepvariant.vcf.gz \\
        --output_gvcf=${sample_id}.deepvariant.g.vcf.gz \\
        --num_shards=${task.cpus} \\
        --intermediate_results_dir ./intermediate_dir
    """
}

process JOINT_CALLING_GLNEXUS {
    tag "GLnexus Joint Calling (${gvcfs.size()} samples)"
    publishDir "${params.outdir}/variants/joint_called", mode: 'copy'

    input:
    path gvcfs

    output:
    path "joint_variants.bcf", emit: bcf
    path "joint_variants.bcf.csi", emit: bcf_index

    script:
    """
    # Create a file with the list of gVCF paths
    ls -1 *.g.vcf.gz > gvcf_list.txt

    # GLnexus is memory intensive, so we set a limit to avoid crashes
    export GLNEXUS_MAX_MEM_MB=\$((${task.memory.toMega()} - 2000))

    glnexus_cli \\
        --config DeepVariantWGS \\
        --list gvcf_list.txt \\
        --dir ./tmp_glnexus \\
    | bcftools view -O b -o joint_variants.bcf -

    bcftools index joint_variants.bcf
    """
}


process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'

_
    input:
    path '*'

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}


/*
 * ====================================================================================
 *   Define Workflow
 * ====================================================================================
 */

workflow {

    // --- Create input channels ---
    Channel
        .fromPath(params.input)
        .splitCsv(header:true, sep:',')
        .map { row -> tuple(row.sample, [ file(row.fastq_1), file(row.fastq_2) ]) }
        .set { ch_fastq_pairs }

    // --- Run pipeline steps ---
    FASTP(ch_fastq_pairs)

    ALIGN_BWA_MEM2(
        FASTP.out.trimmed_reads,
        file(params.ref_fasta),
        file(params.bwa_index)
    )

    MARK_DUPLICATES(ALIGN_BWA_MEM2.out.bam_bai)

    DEEPVARIANT(
        MARK_DUPLICATES.out.dedup_bam_bai,
        file(params.ref_fasta)
    )

    // --- Optional Joint Calling Step ---
    if (params.joint_calling) {
        // Collect all the gVCF files from the DeepVariant output
        // The output tuple is (sample_id, vcf, gvcf, tbi, html). We need the 3rd element (index 2).
        DEEPVARIANT.out.results
            .map { it[2] } // Selects the gvcf path
            .collect()
            .set { ch_all_gvcfs }

        JOINT_CALLING_GLNEXUS(ch_all_gvcfs)
    }


    // --- Collate and run MultiQC ---
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json_report.collect())
    ch_multiqc_files = ch_multiqc_files.mix(ALIGN_BWA_MEM2.out.flagstat.collect())
    // Add more reports to multiqc if available from other tools

    MULTIQC(ch_multiqc_files.collect())
}

workflow.onComplete {
    log.info "Pipeline completed successfully!"
}
