# In this script, we are generating the config file for processing the dnaseq data for heart
import json
heart_ids = []
with open("c://Users/tuyen/Documents/postdoc_asu/projects/PlacentaSexDiff/B_gtex/01_download_data/rna_samples_females_final_forcode.csv", "r") as f:
    for line in f:
        if line.startswith("HeartLeftVentricle"):
            items = line.rstrip("\n").split(",")
            for i in items[1:]:
                j = i.split("-")
                heart_ids.append(j[0] + "-" + j[1])

new_data = {}
new_data["dna_read_group_identifier"] = []
new_data["all_dna_samples"] = []
new_data["all_dna_samples_without_SM"] = []
new_data["dna_samples"] = {}

with open("c://Users/tuyen/Documents/postdoc_asu/projects/PlacentaSexDiff/B_gtex/05_process_dna/process_dna_config.json") as json_file:
    data = json.load(json_file)
    for heart_id in heart_ids:
        for i in data["dna_read_group_identifier"]:
            if i.startswith(heart_id):
                new_data["dna_read_group_identifier"].append(i)
        for i in data:
            if i.startswith(heart_id):
                new_data[i] = data[i]
        for i in data["all_dna_samples"]:
            if i.startswith(heart_id):
                new_data["all_dna_samples"].append(i)
        for i in data["all_dna_samples_without_SM"]:
            if i.startswith(heart_id):
                new_data["all_dna_samples_without_SM"].append(i)
        for i in data["dna_samples"]:
            if i.startswith(heart_id):
                new_data["dna_samples"][i] = data["dna_samples"][i]

new_data["trimmed_fastqs_dna_directory"] = "/scratch/tphung3/PlacentaSexDiff/B_gtex/05_process_dna/trimmed_fastqs_dna/"

with open('c://Users/tuyen/Documents/postdoc_asu/projects/XCI_Methods/B_Heart/01_data_processing/dna_data_processing_config.json', 'w') as outfile:
    json.dump(new_data, outfile)





