configfile: "config.yaml"

# Define a wildcards for sample names if needed, for now, assume single sample
# wildcards:
#     sample="HG002"

# Rule for FASTP
rule fastp:
    input:
        r1=config["input"]["r1"],
        r2=config["input"]["r2"]
    output:
        trimmed_r1=config["outdir"] + "/fastp/HG002_trimmed_R1.fastq.gz",
        trimmed_r2=config["outdir"] + "/fastp/HG002_trimmed_R2.fastq.gz",
        fastp_json=config["outdir"] + "/fastp/HG002_fastp.json",
        fastp_html=config["outdir"] + "/fastp/HG002_fastp.html"
    container:
        config["containers"]["fastp"]
    shell:
        """
        fastp -i {input.r1} -o {output.trimmed_r1} \
              -I {input.r2} -O {output.trimmed_r2} \
              -j {output.fastp_json} -h {output.fastp_html}
        """

# Define the 'all' rule to specify the final output files
# This will be updated as more rules are added
rule all:
    input:
        config["outdir"] + "/fastp/HG002_trimmed_R1.fastq.gz",
        config["outdir"] + "/fastp/HG002_trimmed_R2.fastq.gz",
        config["outdir"] + "/bwa_mem2/HG002.bam",
        config["outdir"] + "/mark_duplicates/HG002.marked.bam",
        config["outdir"] + "/deepvariant/HG002.vcf.gz",
        config["outdir"] + "/deepvariant/HG002.g.vcf.gz",
        # Conditionally include joint calling output
        expand(config["outdir"] + "/glnexus/joint_called.vcf.gz",
               condition=config["joint_calling"]),
        config["outdir"] + "/multiqc/multiqc_report.html"

rule multiqc:
    input:
        fastp_json=rules.fastp.output.fastp_json,
        metrics=rules.mark_duplicates.output.metrics
    output:
        html=config["outdir"] + "/multiqc/multiqc_report.html"
    container:
        config["containers"]["multiqc"]
    resources:
        cpus=1,
        mem_mb=4 * 1024
    shell:
        """
        multiqc -o {config[outdir]}/multiqc {input.fastp_json} {input.metrics}
        """

rule joint_calling_glnexus:
    input:
        gvcf=rules.deepvariant.output.gvcf,
        ref_fasta=config["ref_fasta"]
    output:
        joint_vcf=config["outdir"] + "/glnexus/joint_called.vcf.gz",
        joint_vcf_tbi=config["outdir"] + "/glnexus/joint_called.vcf.gz.tbi"
    container:
        config["containers"]["glnexus"]
    resources:
        cpus=8,
        mem_mb=96 * 1024
    shell:
        """
        glnexus_cli \
            --config DeepVariant \
            --output {output.joint_vcf} \
            {input.gvcf}

        bcftools index -t {output.joint_vcf}
        """

rule deepvariant:
    input:
        marked_bam=rules.mark_duplicates.output.marked_bam,
        ref_fasta=config["ref_fasta"]
    output:
        vcf=config["outdir"] + "/deepvariant/HG002.vcf.gz",
        gvcf=config["outdir"] + "/deepvariant/HG002.g.vcf.gz"
    container:
        config["containers"]["deepvariant"]
    resources:
        cpus=16,
        mem_mb=64 * 1024
    shell:
        """
        run_deepvariant \
            --model_type WGS \
            --ref {input.ref_fasta} \
            --reads {input.marked_bam} \
            --output_vcf {output.vcf} \
            --output_gvcf {output.gvcf} \
            --num_shards {resources.cpus}
        """

rule mark_duplicates:
    input:
        bam=rules.align_bwa_mem2.output.bam
    output:
        marked_bam=config["outdir"] + "/mark_duplicates/HG002.marked.bam",
        metrics=config["outdir"] + "/mark_duplicates/HG002.metrics.txt"
    container:
        config["containers"]["samtools"]
    resources:
        cpus=8,
        mem_mb=32 * 1024
    shell:
        """
        samtools markdup -r -s {input.bam} {output.marked_bam} 2> {output.metrics}
        """

rule align_bwa_mem2:
    input:
        trimmed_r1=rules.fastp.output.trimmed_r1,
        trimmed_r2=rules.fastp.output.trimmed_r2,
        ref_fasta=config["ref_fasta"]
    output:
        bam=config["outdir"] + "/bwa_mem2/HG002.bam"
    container:
        config["containers"]["bwa_mem2"]
    resources:
        cpus=16,
        mem_mb=64 * 1024
    shell:
        """
        bwa-mem2 mem -t {resources.cpus} \
            -K 10000000 \
            -R \"@RG\tID:HG002\tSM:HG002\tLB:HG002\tPL:ILLUMINA\" \
            {input.ref_fasta} \
            {input.trimmed_r1} {input.trimmed_r2} | \
            samtools view -bS - > {output.bam}
        """
