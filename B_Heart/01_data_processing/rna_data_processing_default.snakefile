import os

configfile: "rna_data_processing_config.json"

rule all:
	input:
		expand("processed_bams/rna/{sample}.GRCh38.p12.genome.sorted.bam.bai", sample=config["rna_read_group_identifier"]), #read mapping
		expand("bam_stats/rna/{sample}.GRCh38.p12.genome.rna.sorted.stats", sample=config["rna_read_group_identifier"]), #read mapping
		expand("processed_bams/rna/{sample}.GRCh38.p12.genome.sorted.merged.bam", sample=config["rna_samples"]), #merge
		expand("processed_bams/rna/{sample}.GRCh38.p12.genome.sorted.merged.bam.bai", sample=config["rna_samples"]) #merge

rule hisat2_default:
	input:
		"/data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome/GRCh38.p12.genome.fa"
	output:
		expand(
			"/data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome/GRCh38.p12.genome.{suffix}.ht2",
			suffix=[
				"1", "2", "3", "4", "5", "6", "7", "8"])
	shell:
		"hisat2-build {input} /data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome/GRCh38.p12.genome"

# Read mapping rules
rule hisat2_map_reads:
	input:
		idx = expand(
			"/data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome/GRCh38.p12.genome.{suffix}.ht2",
			suffix=["1", "2", "3", "4", "5", "6", "7", "8"]),
		fq_1 = "trimmed_fastqs_rna/{sample}_trimmed_read1.fastq.gz",
		fq_2 = "trimmed_fastqs_rna/{sample}_trimmed_read2.fastq.gz"
	output:
		"processed_bams/rna/{sample}.GRCh38.p12.genome.sorted.bam"
	threads:
		4
	params:
		threads = 4,
		id = lambda wildcards: config[wildcards.sample]["ID"],
		sm = lambda wildcards: config[wildcards.sample]["SM"],
		lb = lambda wildcards: config[wildcards.sample]["LB"],
		pu = lambda wildcards: config[wildcards.sample]["PU"],
		pl = lambda wildcards: config[wildcards.sample]["PL"]
	shell:
		"PERL5LIB=/home/tphung3/softwares/miniconda3/envs/epitopepipeline/lib/site_perl/5.26.2/ hisat2 -p {params.threads} --dta "
		"--rg-id {params.id} --rg SM:{params.sm} --rg LB:{params.lb} "
		"--rg PU:{params.pu} --rg PL:{params.pl} "
		"-x /data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome/GRCh38.p12.genome "
		"-1 {input.fq_1} -2 {input.fq_2} | "
		"samtools sort -O bam -o {output}"

rule index_bam_rna:
	input:
		"processed_bams/rna/{sample}.GRCh38.p12.genome.sorted.bam"
	output:
		"processed_bams/rna/{sample}.GRCh38.p12.genome.sorted.bam.bai"
	shell:
		"samtools index {input}"

rule bam_stats_rna:
	input:
		"processed_bams/rna/{sample}.GRCh38.p12.genome.sorted.bam"
	output:
		"bam_stats/rna/{sample}.GRCh38.p12.genome.rna.sorted.stats"
	shell:
		"samtools stats {input} | grep ^SN | cut -f 2- > {output}"

rule merge_bams_rna:
	input:
		bams = lambda wildcards: expand(
			"processed_bams/rna/{sample}.GRCh38.p12.genome.sorted.bam",
			sample=config["rna_samples"][wildcards.sample]),
		bais = lambda wildcards: expand(
			"processed_bams/rna/{sample}.GRCh38.p12.genome.sorted.bam.bai",
			sample=config["rna_samples"][wildcards.sample])
	output:
		"processed_bams/rna/{sample}.GRCh38.p12.genome.sorted.merged.bam"
	threads: 4
	params:
		threads = 4
	shell:
		"sambamba merge -t {params.threads} {output} {input.bams}"

rule index_merged_bams_rna:
    input:
        "processed_bams/rna/{sample}.GRCh38.p12.genome.sorted.merged.bam"
    output:
        "processed_bams/rna/{sample}.GRCh38.p12.genome.sorted.merged.bam.bai"
    shell:
        "samtools index {input}"
