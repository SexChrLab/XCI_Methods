import os

configfile: "asereadcounter_config.json"

chromosomes = ["X"]

gatk_path = "/home/tphung3/softwares/gatk-4.1.8.1/gatk"

combiC = []
for key in config["dna_rna"]:
    for item in config["dna_rna"][key]:
        combiC.append((key, item))

import itertools
combiList=list()
for c in combiC:
    combiList.append(c[0]+"_"+c[1])

rule all:
    input:
        expand("asereadcounter/scc/chr{chr}/{combo}_chr{chr}_scc.tsv", combo=combiList, chr=chromosomes)
    input: #Run asereadcounter for autosomes and chrX
        expand("asereadcounter/default/chr{chr}/{combo}_chr{chr}_default.tsv", combo=combiList, chr=chromosomes)

# rule gatk_asereadcounter_default:
#     input:
#         ref = "/data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome/GRCh38.p12.genome.fa",
#         bam = "/scratch/tphung3/XCI_Methods/B_Heart/01_data_processing/processed_bams/rna/{rna}.GRCh38.p12.genome.sorted.merged.bam",
#         sites = "/scratch/tphung3/XCI_Methods/B_Heart/01_data_processing/vqsr/GRCh38.p12.genome.chrX.gatk.called.vqsr.sv.biallelic.snp.{dna}.het.vcf"
#     output:
#         "asereadcounter/default/chr{chr}/{dna}_{rna}_chr{chr}_default.tsv"
#     params:
#         gatk = gatk_path
#     shell:
#         """{params.gatk} ASEReadCounter """
#         """-R {input.ref} """
#         """--output {output} """
#         """--input {input.bam} """
#         """--variant {input.sites} """
#         """--min-depth-of-non-filtered-base 1 """
#         """--min-mapping-quality 10 """
#         """--min-base-quality 10 """

rule gatk_asereadcounter_scc:
    input:
        ref = "/data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome.XXonly/GRCh38.p12.genome.XXonly.fa",
        bam = "/scratch/tphung3/PlacentaSexDiff/B_gtex/04_process_rna/processed_bams/rna/{rna}.GRCh38.p12.genome.XXonly.sorted.merged.bam",
        sites = "/scratch/tphung3/XCI_Methods/B_Heart/01_data_processing/vqsr/GRCh38.p12.genome.XXonly.chrX.gatk.called.vqsr.sv.biallelic.snp.{dna}.het.vcf"
    output:
        "asereadcounter/scc/chr{chr}/{dna}_{rna}_chr{chr}_scc.tsv"
    params:
        gatk = gatk_path
    shell:
        """{params.gatk} ASEReadCounter """
        """-R {input.ref} """
        """--output {output} """
        """--input {input.bam} """
        """--variant {input.sites} """
        """--min-depth-of-non-filtered-base 1 """
        """--min-mapping-quality 10 """
        """--min-base-quality 10 """
