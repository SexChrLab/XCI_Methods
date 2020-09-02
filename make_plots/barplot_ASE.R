# Outside bars
ggplot(data=df, aes(x=dose, y=len)) +
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=len), vjust=-0.3, size=3.5)+
  theme_minimal()
# Inside bars
ggplot(data=df, aes(x=dose, y=len)) +
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=len), vjust=1.6, color="white", size=3.5)+
  theme_minimal()

# Stacked barplot with multiple groups
ggplot(data=df2, aes(x=dose, y=len, fill=supp)) +
  geom_bar(stat="identity")
# Use position=position_dodge()
ggplot(data=df2, aes(x=dose, y=len, fill=supp)) +
  geom_bar(stat="identity", position=position_dodge())

# Change the colors manually
p <- ggplot(data=df2, aes(x=dose, y=len, fill=supp)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()
# Use custom colors
p + scale_fill_manual(values=c('#999999','#E69F00'))
# Use brewer color palettes
p + scale_fill_brewer(palette="Blues")

# read in table
expressedVariants <- read.delim("expressedVariants.txt", sep = "\t", header = TRUE)
expressedVariants$default <- expressedVariants$Overlap + expressedVariants$Unique_default
expressedVariants$Y_masked <- expressedVariants$Overlap + expressedVariants$Unique_SCC
expressedVariants$shared <- expressedVariants$Overlap

library(reshape)
df <- melt(data = expressedVariants, id.vars = c("Sample", "default", "shared", "Y_masked"), measure.vars = c("default", "shared", "Y_masked"))
df$variable <- factor(df$variable, levels = c('B', 'C', 'A'))
library(ggplot2)

q <- ggplot(data=df, aes(x=Sample, y=value, fill=variable)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=value), vjust=1.6, color="white",
            position = position_dodge(0.9), size=3.5)+
 # scale_fill_brewer()+
  scale_fill_manual(values=c("steelblue3", "gray40", "slateblue"))+
  theme_minimal()+
 ggtitle("") +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12)) +
  theme(axis.title=element_text(size=15)) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=15)) +
  labs(y = "Expressed heterozygous sites")+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1.25, hjust=1, margin=margin(0,0,0,0)))
q
library("gridExtra")
ggsave("expressedHetSites_barplot.pdf", grid.arrange(q), height = 5, width = 8)

          