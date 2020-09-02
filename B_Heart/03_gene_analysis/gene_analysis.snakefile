import os

configfile: "gene_analysis_config.json"

rule all:
    input:
        expand("results/{mapping_type}/chrX/{sample}_chrX_{mapping_type}_geneinfo_unique.tsv", sample=config["all_samples"], mapping_type=config["mapping_types"])

rule convert_asereadcounter_to_bed:
    input:
        "/scratch/tphung3/XCI_Methods/B_Heart/02_run_asereadcounter/asereadcounter/{mapping_type}/chrX/{sample}_chrX_{mapping_type}.tsv"
    output:
        "results/{mapping_type}/chrX/{sample}_chrX_{mapping_type}.bed"
    shell:
        """
        python /scratch/tphung3/PlacentaSexDiff/E_escape_genes/wes_genotyping/convert_asereadcounter_to_bed.py --input_asereadcounter {input} --output_bed {output}
        """
rule bedtools_intersect:
    input:
        "results/{mapping_type}/chrX/{sample}_chrX_{mapping_type}.bed"
    output:
        "results/{mapping_type}/chrX/{sample}_chrX_{mapping_type}_geneinfo.tsv"
    shell:
        """
        bedtools intersect -a /scratch/tphung3/PlacentaSexDiff/E_escape_genes/wes_genotyping/gtf_bed/gencode.v29.annotation.chrX.bed -b {input} -wa -wb > {output}
        """

rule find_unique_lines_after_bedtools:
    input:
        "results/{mapping_type}/chrX/{sample}_chrX_{mapping_type}_geneinfo.tsv"
    output:
        "results/{mapping_type}/chrX/{sample}_chrX_{mapping_type}_geneinfo_unique.tsv"
    shell:
        """
        python /scratch/tphung3/PlacentaSexDiff/E_escape_genes/wes_genotyping/find_unique_lines_after_bedtools.py --infile {input} --outfile {output} --index 5
        """
