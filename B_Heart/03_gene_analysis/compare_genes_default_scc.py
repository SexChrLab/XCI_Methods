import sys

default_genes = set()
with open(sys.argv[1], "r") as f:
    for line in f:
        items = line.rstrip("\n").split("\t")
        if int(items[12]) > 10:
            default_genes.add(items[3])

scc_genes = set()
with open(sys.argv[2], "r") as f:
    for line in f:
        items = line.rstrip("\n").split("\t")
        if int(items[12]) > 10:
            scc_genes.add(items[3])

# print(len(default_genes.intersection(scc_genes)))
# print(len(default_genes-scc_genes))
# print(len(scc_genes-default_genes))

print(default_genes-scc_genes)
print(scc_genes-default_genes)
