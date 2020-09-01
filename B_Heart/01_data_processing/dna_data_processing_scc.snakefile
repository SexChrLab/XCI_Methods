import os

configfile: "dna_data_processing_config.json"

sex_chr = ["X"]
gatk_path = "/home/tphung3/softwares/gatk-4.1.8.1/gatk"

rule all:
    input:
        "vqsr/GRCh38.p12.genome.XXonly.chrX.gatk.called.vqsr.sv.vcf.gz",
        expand("vqsr/GRCh38.p12.genome.XXonly.chrX.gatk.called.vqsr.sv.biallelic.snp.{sample}.vcf", sample=config["all_dna_samples_without_SM"]),
        expand("vqsr/GRCh38.p12.genome.XXonly.chrX.gatk.called.vqsr.sv.biallelic.snp.{sample}.het.vcf", sample=config["all_dna_samples_without_SM"])


rule gatk_combinegvcfs_chrX:
	input:
		ref = "/data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome.XXonly/GRCh38.p12.genome.XXonly.fa",
		gvcfs = expand(
			"/scratch/tphung3/PlacentaSexDiff/B_gtex/05_process_dna/gvcfs/gvcfs/{sample}/{sample}.GRCh38.p12.genome.XXonly.chrX.g.vcf.gz", sample=config["all_dna_samples"])
	output:
		"combined_gvcfs/GRCh38.p12.genome.XXonly.chrX.gatk.combinegvcf.g.vcf.gz"
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
		ref = "/data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome.XXonly/GRCh38.p12.genome.XXonly.fa",
		gvcf = "combined_gvcfs/GRCh38.p12.genome.XXonly.chrX.gatk.combinegvcf.g.vcf.gz"
	output:
		"genotyped_vcfs/GRCh38.p12.genome.XXonly.chrX.gatk.called.raw.vcf.gz"
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
        ref = "/data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome.XXonly/GRCh38.p12.genome.XXonly.fa",
        vcf = "genotyped_vcfs/GRCh38.p12.genome.XXonly.chrX.gatk.called.raw.vcf.gz",
        hapmap = "/data/CEM/shared/public_data/validated_variant_resources/hapmap_3.3.hg38.vcf.gz",
        omni = "/data/CEM/shared/public_data/validated_variant_resources/1000G_omni2.5.hg38.vcf.gz",
        thousandG = "/data/CEM/shared/public_data/validated_variant_resources/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
    output:
        recal = "vqsr/GRCh38.p12.genome.XXonly.chrX_output.recal",
        tranches = "vqsr/GRCh38.p12.genome.XXonly.chrX_output.tranches"
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
        ref = "/data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome.XXonly/GRCh38.p12.genome.XXonly.fa",
        vcf = "genotyped_vcfs/GRCh38.p12.genome.XXonly.chrX.gatk.called.raw.vcf.gz",
        tranches = "vqsr/GRCh38.p12.genome.XXonly.chrX_output.tranches",
        recal = "vqsr/GRCh38.p12.genome.XXonly.chrX_output.recal"
    output:
        "vqsr/GRCh38.p12.genome.XXonly.chrX.gatk.called.vqsr.vcf.gz"
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
        ref = "/data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome.XXonly/GRCh38.p12.genome.XXonly.fa",
        vcf = "vqsr/GRCh38.p12.genome.XXonly.chrX.gatk.called.vqsr.vcf.gz"
    output:
        "vqsr/GRCh38.p12.genome.XXonly.chrX.gatk.called.vqsr.sv.vcf.gz"
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
        ref = "/data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome.XXonly/GRCh38.p12.genome.XXonly.fa",
        vcf = "vqsr/GRCh38.p12.genome.XXonly.chrX.gatk.called.vqsr.sv.vcf.gz"
    output:
        "vqsr/GRCh38.p12.genome.XXonly.chrX.gatk.called.vqsr.sv.biallelic.snp.vcf.gz"
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
#         "vqsr/GRCh38.p12.genome.XXonly.chrX.gatk.called.vqsr.sv.biallelic.snp.vcf.gz"
#     output:
#         "vqsr/GRCh38.p12.genome.XXonly.chrX.gatk.called.vqsr.sv.biallelic.snp.{sample}.vcf"
#     params:
#         sample = "{sample}"
#     shell:
#         """bcftools view -s {params.sample} {input} > {output}"""

# After subsetting for each individual. In some individuals,
# the genotypes could be homozygous for the reference. This next rule is to remove these sites.
rule gatk_selectheterozygous:
    input:
        ref = "/data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome.XXonly/GRCh38.p12.genome.XXonly.fa",
        vcf = "vqsr/GRCh38.p12.genome.XXonly.chrX.gatk.called.vqsr.sv.biallelic.snp.{sample}.vcf"
    output:
        "vqsr/GRCh38.p12.genome.XXonly.chrX.gatk.called.vqsr.sv.biallelic.snp.{sample}.het.vcf"
    params:
        gatk = gatk_path
    shell:
        """{params.gatk} SelectVariants """
        """-R {input.ref} """
        """-V {input.vcf} """
        """-O {output} """
        """-select "AC == 1" """
