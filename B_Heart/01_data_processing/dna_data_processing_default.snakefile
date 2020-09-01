import os

configfile: "dna_data_processing_config.json"

sex_chr = ["X"]
gatk_path = "/home/tphung3/softwares/gatk-4.1.8.1/gatk"

# ------------------------
# Default reference genome
# ------------------------
rule all:
    input:
        "/data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome/GRCh38.p12.genome.fa.fai",
        "/data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome/GRCh38.p12.genome.dict",
        "/data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome/GRCh38.p12.genome.fa.amb",
        expand("processed_bams/dna/default/{sample}.GRCh38.p12.genome.sorted.bam.bai", sample=config["dna_read_group_identifier"]),
        expand("bam_stats/dna/default/{sample}.GRCh38.p12.genome.sorted.stats", sample=config["dna_read_group_identifier"]),
        expand("processed_bams/dna/default/{sample}.GRCh38.p12.genome.sorted.merged.bam", sample=config["all_dna_samples"]),
        expand("processed_bams/dna/default/{sample}.GRCh38.p12.genome.sorted.merged.bam.bai", sample=config["all_dna_samples"]),
        expand("picard_stats/dna/default/{sample}.GRCh38.p12.genome.picard_mkdup_metrics.txt", sample=config["all_dna_samples"]),
        expand("processed_bams/dna/default/{sample}.GRCh38.p12.genome.sorted.merged.mkdup.bam.bai", sample=config["all_dna_samples"]),
        expand("gvcfs/{sample}/{sample}.GRCh38.p12.genome.chr{chr}.g.vcf.gz", sample=config["all_dna_samples"], chr=sex_chr),
        "combined_gvcfs/GRCh38.p12.genome.chrX.gatk.combinegvcf.g.vcf.gz",
        expand("genotyped_vcfs/GRCh38.p12.genome.chr{chr}.gatk.called.raw.vcf.gz", chr=sex_chr),
        "vqsr/GRCh38.p12.genome.chrX.gatk.called.vqsr.sv.vcf.gz",
        expand("vqsr/GRCh38.p12.genome.chrX.gatk.called.vqsr.sv.biallelic.snp.{sample}.vcf", sample=config["all_dna_samples_without_SM"]),
        expand("vqsr/GRCh38.p12.genome.chrX.gatk.called.vqsr.sv.biallelic.snp.{sample}.het.vcf", sample=config["all_dna_samples_without_SM"])



rule prep_default_ref:
	input:
		"/data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome/GRCh38.p12.genome.fa"
	output:
		fai = "/data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome/GRCh38.p12.genome.fa.fai",
		dict = "/data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome/GRCh38.p12.genome.dict",
		amb = "/data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome/GRCh38.p12.genome.fa.amb"
	shell:
		"""
		samtools faidx {input};
		samtools dict -o {output.dict} {input};
		bwa index {input}
		"""

rule map_default:
	input:
		fq_1 = os.path.join(config["trimmed_fastqs_dna_directory"], "{sample}_trimmed_read1.fastq.gz"),
		fq_2 = os.path.join(config["trimmed_fastqs_dna_directory"], "{sample}_trimmed_read2.fastq.gz"),
		fai = "/data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome/GRCh38.p12.genome.fa.fai",
		ref = "/data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome/GRCh38.p12.genome.fa"
	output:
		"processed_bams/dna/default/{sample}.GRCh38.p12.genome.sorted.bam"
	params:
		id = lambda wildcards: config[wildcards.sample]["ID"],
		sm = lambda wildcards: config[wildcards.sample]["SM"],
		lb = lambda wildcards: config[wildcards.sample]["LB"],
		pu = lambda wildcards: config[wildcards.sample]["PU"],
		pl = lambda wildcards: config[wildcards.sample]["PL"]
	threads:
		4
	shell:
		"bwa mem -t {threads} -R "
		"'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
		"{input.ref} {input.fq_1} {input.fq_2} "
		"| samblaster "
		"| samtools fixmate -O bam - - | samtools sort "
		"-O bam -o {output}"

rule index_bam_dna_default:
	input:
		"processed_bams/dna/default/{sample}.GRCh38.p12.genome.sorted.bam"
	output:
		"processed_bams/dna/default/{sample}.GRCh38.p12.genome.sorted.bam.bai"
	shell:
		"samtools index {input}"

rule bam_stats_dna_default:
	input:
		"processed_bams/dna/default/{sample}.GRCh38.p12.genome.sorted.bam"
	output:
		"bam_stats/dna/default/{sample}.GRCh38.p12.genome.sorted.stats"
	shell:
		"samtools stats {input} | grep ^SN | cut -f 2- > {output}"

rule merge_bams_dna_default:
	input:
		bams = lambda wildcards: expand(
			"processed_bams/dna/default/{sample}.GRCh38.p12.genome.sorted.bam",
			sample=config["dna_samples"][wildcards.sample]),
		bais = lambda wildcards: expand(
			"processed_bams/dna/default/{sample}.GRCh38.p12.genome.sorted.bam.bai",
			sample=config["dna_samples"][wildcards.sample])
	output:
		"processed_bams/dna/default/{sample}.GRCh38.p12.genome.sorted.merged.bam"
	threads: 4
	params:
		threads = 4
	shell:
		"sambamba merge -t {params.threads} {output} {input.bams}"

rule index_merged_bams_dna_default:
    input:
        "processed_bams/dna/default/{sample}.GRCh38.p12.genome.sorted.merged.bam"
    output:
        "processed_bams/dna/default/{sample}.GRCh38.p12.genome.sorted.merged.bam.bai"
    shell:
        "samtools index {input}"

rule picard_mkdups_default:
    input:
        bam = "processed_bams/dna/default/{sample}.GRCh38.p12.genome.sorted.merged.bam",
        bai = "processed_bams/dna/default/{sample}.GRCh38.p12.genome.sorted.merged.bam.bai"
    output:
        bam = "processed_bams/dna/default/{sample}.GRCh38.p12.genome.sorted.merged.mkdup.bam",
        metrics = "picard_stats/dna/default/{sample}.GRCh38.p12.genome.picard_mkdup_metrics.txt"
    threads: 4
    shell:
        "picard -Xmx14g MarkDuplicates I={input.bam} O={output.bam} "
        "M={output.metrics}"

rule index_mkdup_bams_dna_default:
    input:
        "processed_bams/dna/default/{sample}.GRCh38.p12.genome.sorted.merged.mkdup.bam"
    output:
        "processed_bams/dna/default/{sample}.GRCh38.p12.genome.sorted.merged.mkdup.bam.bai"
    shell:
        "samtools index {input}"

# ----------------------------------
# Joint-genotype on the X chromosome
# ----------------------------------
rule gatk_gvcf:
	input:
		ref = "/data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome/GRCh38.p12.genome.fa",
		bam = "processed_bams/dna/default/{sample}.GRCh38.p12.genome.sorted.merged.mkdup.bam",
		bai = "processed_bams/dna/default/{sample}.GRCh38.p12.genome.sorted.merged.mkdup.bam.bai"
	output:
		"gvcfs/{sample}/{sample}.GRCh38.p12.genome.chr{chr}.g.vcf.gz"
	params:
		gatk = gatk_path,
		chrm_n = "chr{chr}"
	threads:
		4
	shell:
		"{params.gatk} "
		"HaplotypeCaller -R {input.ref} -I {input.bam} -L {params.chrm_n} "
		"--emit-ref-confidence GVCF --output {output}"

rule gatk_combinegvcfs_chrX:
	input:
		ref = "/data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome/GRCh38.p12.genome.fa",
		gvcfs = expand(
			"gvcfs/{sample}/{sample}.GRCh38.p12.genome.chrX.g.vcf.gz", sample=config["all_dna_samples"])
	output:
		"combined_gvcfs/GRCh38.p12.genome.chrX.gatk.combinegvcf.g.vcf.gz"
	params:
		gatk = gatk_path,
		chr_n = "chrX"
	threads:
		4
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""{params.gatk} --java-options "-Xmx10g" """
			"""CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chr_n} -O {output}""")

rule gatk_genotypegvcf:
	input:
		ref = "/data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome/GRCh38.p12.genome.fa",
		gvcf = "combined_gvcfs/GRCh38.p12.genome.chr{chr}.gatk.combinegvcf.g.vcf.gz"
	output:
		"genotyped_vcfs/GRCh38.p12.genome.chr{chr}.gatk.called.raw.vcf.gz"
	params:
		gatk = gatk_path
	threads:
		4
	shell:
		"""{params.gatk} --java-options "-Xmx10g" """
		"""GenotypeGVCFs -R {input.ref} -V {input.gvcf} -O {output}"""

# chrX
rule gatk_variantrecalibrator_chrX:
    input:
        ref = "/data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome/GRCh38.p12.genome.fa",
        vcf = "genotyped_vcfs/GRCh38.p12.genome.chrX.gatk.called.raw.vcf.gz",
        hapmap = "/data/CEM/shared/public_data/validated_variant_resources/hapmap_3.3.hg38.vcf.gz",
        omni = "/data/CEM/shared/public_data/validated_variant_resources/1000G_omni2.5.hg38.vcf.gz",
        thousandG = "/data/CEM/shared/public_data/validated_variant_resources/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
    output:
        recal = "vqsr/GRCh38.p12.genome.chrX_output.recal",
        tranches = "vqsr/GRCh38.p12.genome.chrX_output.tranches"
    params:
        gatk = gatk_path
    shell:
        """{params.gatk} --java-options "-Xmx16g" VariantRecalibrator """
        """-R {input.ref} -V {input.vcf}  """
        """--resource:hapmap,known=false,training=true,truth=true,prior=15.0 {input.hapmap} """
        """--resource:omni,known=false,training=true,truth=false,prior=12.0 {input.omni} """
        """--resource:1000G,known=false,training=true,truth=false,prior=10.0 {input.thousandG} """
        """-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff """
        """-mode SNP """
        """-O {output.recal} """
        """--tranches-file {output.tranches} """

rule gatk_applyvqsr_chrX:
    input:
        ref = "/data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome/GRCh38.p12.genome.fa",
        vcf = "genotyped_vcfs/GRCh38.p12.genome.chrX.gatk.called.raw.vcf.gz",
        tranches = "vqsr/GRCh38.p12.genome.chrX_output.tranches",
        recal = "vqsr/GRCh38.p12.genome.chrX_output.recal"
    output:
        "vqsr/GRCh38.p12.genome.chrX.gatk.called.vqsr.vcf.gz"
    params:
        gatk = gatk_path
    shell:
        """{params.gatk} --java-options "-Xmx16g" ApplyVQSR """
        """-R {input.ref} """
        """-V {input.vcf} """
        """-O {output} """
        """--truth-sensitivity-filter-level 99.0 """
        """--tranches-file {input.tranches} """
        """--recal-file {input.recal} """
        """-mode SNP """

rule gatk_selectvariants_chrX:
    input:
        ref = "/data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome/GRCh38.p12.genome.fa",
        vcf = "vqsr/GRCh38.p12.genome.chrX.gatk.called.vqsr.vcf.gz"
    output:
        "vqsr/GRCh38.p12.genome.chrX.gatk.called.vqsr.sv.vcf.gz"
    params:
        gatk = gatk_path
    shell:
        """{params.gatk} --java-options "-Xmx16g" SelectVariants """
        """-R {input.ref} """
        """-V {input.vcf} """
        """--exclude-filtered """
        """-O {output} """

#----------------------------------------
# Further processing VCF
# 1. Restrict to biallelic sites
# 2. Subset VCF files for each individual
# 3. Keep only the heterozygous sites
# Do this for chr8 and chrX
#----------------------------------------
rule gatk_selectbiallelic:
    input:
        ref = "/data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome/GRCh38.p12.genome.fa",
        vcf = "vqsr/GRCh38.p12.genome.chrX.gatk.called.vqsr.sv.vcf.gz"
    output:
        "vqsr/GRCh38.p12.genome.chrX.gatk.called.vqsr.sv.biallelic.snp.vcf.gz"
    params:
        gatk = gatk_path
    shell:
        """{params.gatk} SelectVariants """
        """-R {input.ref} """
        """-V {input.vcf} """
        """-O {output} """
        """--select-type-to-include SNP """
        """--restrict-alleles-to BIALLELIC """

# rule subset_individuals_vqsr:
#     input:
#         "vqsr/GRCh38.p12.genome.chrX.gatk.called.vqsr.sv.biallelic.snp.vcf.gz"
#     output:
#         "vqsr/GRCh38.p12.genome.chrX.gatk.called.vqsr.sv.biallelic.snp.{sample}.vcf"
#     params:
#         sample = "{sample}"
#     shell:
#         """bcftools view -s {params.sample} {input} > {output}"""

# After subsetting for each individual. In some individuals,
# the genotypes could be homozygous for the reference. This next rule is to remove these sites.
rule gatk_selectheterozygous:
    input:
        ref = "/data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome/GRCh38.p12.genome.fa",
        vcf = "vqsr/GRCh38.p12.genome.chrX.gatk.called.vqsr.sv.biallelic.snp.{sample}.vcf"
    output:
        "vqsr/GRCh38.p12.genome.chrX.gatk.called.vqsr.sv.biallelic.snp.{sample}.het.vcf"
    params:
        gatk = gatk_path
    shell:
        """{params.gatk} SelectVariants """
        """-R {input.ref} """
        """-V {input.vcf} """
        """-O {output} """
        """-select "AC == 1" """
