library(ggplot2)

# Counter({'RTL8B': 4, 'EMD': 3, 'MAGED1': 3, 'EIF2S3': 2, 'WDR45': 2, 'SMARCA1': 1, 'ITM2A': 1, 'APOO': 1, 'UXT': 1, 'ARMCX5-GPRASP2': 1, 'DNASE1L1': 1, 'AF241726.2': 1, 'FHL1': 1})
# Counter({'DHRSX': 9, 'ASMTL': 9, 'GTPBP6': 8, 'CD99': 7, 'RBMX': 5, 'AKAP17A': 5, 'ASMTL-AS1': 5, 'TIMM17B': 4, 'IL3RA': 3, 'MXRA5': 3, 'RBM3': 2, 'WWC3': 2, 'PCDH11X': 1, 'ARSD': 1, 'IL1RAPL2': 1, 'IGSF1': 1, 'AL034405.1': 1, 'TSIX': 1, 'SLC25A6': 1, 'WASH6P': 1})


# -------
# default
# -------
default = data.frame(genes = c('RTL8B', 'EMD', 'MAGED1', 'EIF2S3', 'WDR45', 'SMARCA1', 'ITM2A', 'APOO', 'UXT', 'ARMCX5-GPRASP2', 'DNASE1L1', 'AF241726.2', 'FHL1'), counts = c(4, 3, 3, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1))

default$genes = as.character(default$genes)
default$genes = factor(default$genes, levels = unique(default$genes))

ggplot(data = default, aes(x=genes, y=counts)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(x="Genes that are unique to default reference", y="Counts") +
  theme(axis.text.x = element_text(angle=25)) +
  coord_cartesian(ylim=c(0,10))
ggsave("c://Users/tuyen/Documents/postdoc_asu/projects/XCI_Methods/B_Heart/03_gene_analysis/plots/default_unique_genes.png", width = 7, height = 3.5, units = "in")

# ---
# scc
# ---
scc = data.frame(genes = c('DHRSX', 'ASMTL', 'GTPBP6', 'CD99', 'RBMX', 'AKAP17A', 'ASMTL-AS1', 'TIMM17B', 'IL3RA', 'MXRA5', 'RBM3', 'WWC3', 'PCDH11X', 'ARSD', 'IL1RAPL2', 'IGSF1', 'AL034405.1', 'TSIX', 'SLC25A6', 'WASH6P'), counts = c(9, 9, 8, 7, 5, 5, 5, 4, 3, 3, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1))

scc$genes = as.character(scc$genes)
scc$genes = factor(scc$genes, levels = unique(scc$genes))

ggplot(data = scc, aes(x=genes, y=counts)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(x="Genes that are unique to scc reference", y="Counts") +
  theme(axis.text.x = element_text(angle=25)) +
  coord_cartesian(ylim=c(0,10))
ggsave("c://Users/tuyen/Documents/postdoc_asu/projects/XCI_Methods/B_Heart/03_gene_analysis/plots/scc_unique_genes.png", width = 7, height = 3.5, units = "in")