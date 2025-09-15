# snakefile

import pandas as pd

# --- Load Configuration ---
# All parameters, paths, and container URIs are defined in the config.yaml file.
configfile: "config.yaml"

# --- Load Sample Sheet ---
# Read the sample information from a CSV file into a pandas DataFrame.
# The sample names are used as wildcards throughout the pipeline.
SAMPLES_DF = pd.read_csv(config["samplesheet"], sep=",").set_index("sample", drop=False)
SAMPLES = list(SAMPLES_DF.index)


# --- Target Rule (rule all) ---
# This is the main entry point of the workflow. Snakemake will determine which
# rules need to be run to generate the files listed here. This process of
# working backward from the desired output is a core principle of Snakemake.
rule all:
    input:
        # Collect MultiQC report which depends on all other analysis steps.
        "results/multiqc/multiqc_report.html",
        # Optionally, include the joint calling results if the flag is set in the config.
        "results/variants/joint_called/joint_variants.bcf" if config["joint_calling"] else []

# ====================================================================================
#   Define Analysis Rules
# ====================================================================================

# --- Helper Function for Input Files ---
# This function retrieves the paths to the FASTQ files for a given sample.
# Using a helper function can make the rule inputs cleaner and more readable.
def get_fastqs(wildcards):
    return {
        'r1': SAMPLES_DF.loc[wildcards.sample, 'fastq_1'],
        'r2': SAMPLES_DF.loc[wildcards.sample, 'fastq_2']
    }

# --- Rule: FASTP ---
# Performs quality control and trimming on raw FASTQ reads.
rule fastp:
    input:
        unpack(get_fastqs)
    output:
        r1 = "results/fastp/{sample}/{sample}.trimmed.R1.fastq.gz",
        r2 = "results/fastp/{sample}/{sample}.trimmed.R2.fastq.gz",
        html = "results/fastp/{sample}/{sample}.fastp.html",
        json = "results/fastp/{sample}/{sample}.fastp.json"
    params:
        # Additional fastp parameters can be specified in the config file.
        extra = config.get("params", {}).get("fastp", "")
    log:
        "logs/fastp/{sample}.log"
    threads: 64
    container:
        config["containers"]["fastp"]
    shell:
        """
        fastp \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.r1} \
            -O {output.r2} \
            -h {output.html} \
            -j {output.json} \
            -w {threads} \
            {params.extra} >& {log}
        """

# --- Rule: ALIGN_BWA_MEM2 ---
# Aligns trimmed reads to the reference genome.
# The reference genome and its index are marked as `ancient` files, meaning
# Snakemake will not attempt to regenerate them if they are older than the input.
rule align_bwa_mem2:
    input:
        r1 = "results/fastp/{sample}/{sample}.trimmed.R1.fastq.gz",
        r2 = "results/fastp/{sample}/{sample}.trimmed.R2.fastq.gz",
        ref_fasta = ancient(config["ref_fasta"]),
        bwa_index = ancient(config["bwa_index"]) # This should be the prefix, not the files themselves
    output:
        bam = "results/bam/{sample}/{sample}.align.bam"
    params:
        # The read group string is essential for downstream tools.
        rg_string = lambda wildcards: f"@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}\tPL:ILLUMINA",
    log:
        "logs/align_bwa_mem2/{sample}.log"
    threads: 64
    container:
        config["containers"]["bwa"]
    shell:
        """
        (bwa-mem2 mem \
            -t {threads} \
            -R '{params.rg_string}' \
            {input.bwa_index} \
            {input.r1} \
            {input.r2} | samtools view -b -o {output.bam}) >& {log}
        """

# --- Rule: BAM_PROCESSING ---
# This rule combines several steps to process the aligned SAM file:
# 1. Sorts the SAM file and converts it to BAM format.
# 2. Marks PCR duplicates to improve variant calling accuracy.
# 3. Indexes the final BAM file for fast access.
# Chaining these commands with pipes is more efficient as it avoids writing
# intermediate files to disk.
rule bam_processing:
    input:
        bam = "results/bam/{sample}/{sample}.align.bam"
    output:
        bam = "results/bam/{sample}/{sample}.dedup.bam",
        bai = "results/bam/{sample}/{sample}.dedup.bam.bai",
        flagstat = "results/bam/{sample}/{sample}.flagstat.txt"
    log:
        "logs/bam_processing/{sample}.log"
    threads: 64
    container:
        config["containers"]["samtools"]
    shell:
        """
        (samtools collate -@ 8 -O {input.bam} | \
            samtools fixmate -m -@ 8 - - | \
            samtools sort -@ 8 - | \
            samtools markdup -r -@ 8 - {output.bam}) >& {log}

        samtools index {output.bam}
        samtools flagstat {input.bam} > {output.flagstat}
        """



# --- Rule: DEEPVARIANT ---
# Calls variants for a single sample using the DeepVariant model.
# It uses a temporary directory for intermediate files, which is automatically
# cleaned up by Snakemake upon successful rule completion.
rule deepvariant:
    input:
        bam = "results/bam/{sample}/{sample}.dedup.bam",
        bai = "results/bam/{sample}/{sample}.dedup.bam.bai",
        ref_fasta = config["ref_fasta"]
    output:
        vcf = "results/variants/single_sample/{sample}/{sample}.deepvariant.vcf.gz",
        gvcf = "results/variants/single_sample/{sample}/{sample}.deepvariant.g.vcf.gz",
        tbi = "results/variants/single_sample/{sample}/{sample}.deepvariant.vcf.gz.tbi",
        html = "results/variants/single_sample/{sample}/{sample}.deepvariant.visual_report.html",
        intermediate_dir = temp("intermediate_dir/{sample}")
    log:
        "logs/deepvariant/{sample}.log"
    threads: 64
    container:
        config["containers"]["deepvariant"]
    shell:
        """
        /opt/deepvariant/bin/run_deepvariant \
            --model_type=WGS \
            --ref=$(realpath {input.ref_fasta}) \
            --reads={input.bam} \
            --output_vcf={output.vcf} \
            --output_gvcf={output.gvcf} \
            --num_shards={threads} \
            --vcf_stats_report=true \
            --intermediate_results_dir {output.intermediate_dir} >& {log}

        # The HTML report is generated inside the intermediate directory.
        # We move it to the final output location.
        mv {output.intermediate_dir}/report.html {output.html}
        """


# --- Rule: JOINT_CALLING_GLNEXUS ---
# Aggregates gVCFs from all samples and performs joint variant calling.
# This is an "aggregation" rule that runs only once.
rule joint_calling_glnexus:
    input:
        # Use expand to gather all g.vcf.gz files from the deepvariant rule.
        gvcfs = expand("results/variants/single_sample/{sample}/{sample}.deepvariant.g.vcf.gz", sample=SAMPLES)
    output:
        bcf = "results/variants/joint_called/joint_variants.bcf",
        csi = "results/variants/joint_called/joint_variants.bcf.csi"
    log:
        "logs/joint_calling/glnexus.log"
    threads: 64
    resources:
        # GLnexus is memory intensive. Define memory as a resource.
        mem_mb=100000
    container:
        config["containers"]["glnexus"]
    shell:
        """
        # Create a file containing the paths to all input gVCFs.
        # Using `echo` is more robust than `ls`.
        for gvcf in {input.gvcfs}; do echo $gvcf; done > gvcf_list.txt

        # Set memory limit for GLnexus to prevent crashes.
        export GLNEXUS_MAX_MEM_MB=$(({resources.mem_mb} - 2000))

        (glnexus_cli \\
            --config DeepVariantWGS \\
            --list gvcf_list.txt \\
            --dir ./tmp_glnexus \\
        | bcftools view -O b -o {output.bcf} -) >& {log}

        bcftools index {output.bcf}
        """

# --- Rule: MULTIQC ---
# Gathers reports from various tools and creates a single summary HTML file.
# This is the final aggregation rule.
rule multiqc:
    input:
        # Use expand to collect all fastp json reports and bwa flagstat files.
        fastp_json = expand("results/fastp/{sample}/{sample}.fastp.json", sample=SAMPLES),
        flagstat = expand("results/bam/{sample}/{sample}.flagstat.txt", sample=SAMPLES)
    output:
        report = "results/multiqc/multiqc_report.html"
    log:
        "logs/multiqc/multiqc.log"
    container:
        config["containers"]["multiqc"]
    shell:
        """
        # MultiQC scans the specified directory for known report files.
        # We'll run it in the top-level 'results' directory to catch everything.
        multiqc results/ --outdir results/multiqc -f >& {log}
        """

