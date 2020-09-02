import json
heart_ids = []
heart_ids_set = set()
with open("c://Users/tuyen/Documents/postdoc_asu/projects/PlacentaSexDiff/B_gtex/01_download_data/rna_samples_females_final_forcode.csv", "r") as f:
    for line in f:
        if line.startswith("HeartLeftVentricle"):
            items = line.rstrip("\n").split(",")
            for i in items[1:]:
                j = i.split("-")
                heart_ids.append(j[0] + "-" + j[1])
                heart_ids_set.add(i)

new_data = {}
new_data["dna_rna"] = {}

with open("c://Users/tuyen/Documents/postdoc_asu/projects/PlacentaSexDiff/B_gtex/06_run_asereadcounter/asereadcounter_config.json") as json_file:
    data = json.load(json_file)
    for heart_id in heart_ids:
        for i in data["dna_rna"]:
            if i.startswith(heart_id):
                for j in data["dna_rna"][i]:
                    if j in heart_ids_set:
                        new_data["dna_rna"][i] = [j]

with open('c://Users/tuyen/Documents/postdoc_asu/projects/XCI_Methods/B_Heart/02_run_asereadcounter/asereadcounter_config.json', 'w') as outfile:
    json.dump(new_data, outfile)