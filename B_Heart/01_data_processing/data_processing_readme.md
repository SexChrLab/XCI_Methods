## Whole exome
1. Map to a default reference genome and to a sex-chromosome-complement (SCC) reference genome
2. Start with the trimmed reads (TODO: add more information here)
3. Generate the config file (`dna_data_processing_config.json`) by running the Python script `generate_json_config_dna.py`
4. Default genome: `data_dna_processing_default.snakefile`
5. SCC genome: `data_dna_processing_scc.snakefile`

## RNAseq
1. Map to a default reference genome
2. Generate the config file (`rna_data_processing_config.json`) by running the Python script `generate_json_config_rna.py`
3. Default genome: `rna_data_processing_default.snakefile`
