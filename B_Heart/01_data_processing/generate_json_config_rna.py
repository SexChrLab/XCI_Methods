# In this script, we are generating the config file for processing the rnaseq data for heart
import json
heart_ids = []
with open("c://Users/tuyen/Documents/postdoc_asu/projects/PlacentaSexDiff/B_gtex/01_download_data/rna_samples_females_final_forcode.csv", "r") as f:
    for line in f:
        if line.startswith("HeartLeftVentricle"):
            items = line.rstrip("\n").split(",")
            for i in items[1:]:
                heart_ids.append(i)

new_data = {}
new_data["rna_read_group_identifier"] = []
new_data["all_rna_samples"] = []
new_data["rna_samples"] = {}

with open("c://Users/tuyen/Documents/postdoc_asu/projects/PlacentaSexDiff/B_gtex/04_process_rna/process_rna_config.json") as json_file:
    data = json.load(json_file)
    for heart_id in heart_ids:
        for i in data["rna_read_group_identifier"]:
            if i.startswith(heart_id):
                new_data["rna_read_group_identifier"].append(i)
        for i in data:
            if i.startswith(heart_id):
                new_data[i] = data[i]
        for i in data["all_rna_samples"]:
            if i.startswith(heart_id):
                new_data["all_rna_samples"].append(i)
        for i in data["rna_samples"]:
            if i.startswith(heart_id):
                new_data["rna_samples"][i] = data["rna_samples"][i]

with open('c://Users/tuyen/Documents/postdoc_asu/projects/XCI_Methods/B_Heart/01_data_processing/rna_data_processing_config.json', 'w') as outfile:
    json.dump(new_data, outfile)





