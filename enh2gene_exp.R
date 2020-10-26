###### are enhancers in dv groups functional?

## very general level...
library(regioneR)
library(ggplot2)

## create list of CRE groups
dfs.gene.exp <- list(df7.all, df8.all, df9)

## prepare RNA-seq data with coordinates extended up and downstream of each TSS
tss.ext.coord <- tss.ext.qlf.cut[, 1:4]
colnames(tss.ext.coord) <- c("chr", "start", "end", "logFC")
tss.ext.qlf.gr <- toGRanges(tss.ext.coord)

## iterate through CRE groups and get log2FC for each transcript overlapping CREs
gene.exp.dv <- data.frame()

for(i in dfs.gene.exp){
  df.gr <- toGRanges(i)
  ols <- findOverlaps(df.gr, tss.ext.qlf.gr)
  gene.exp <- tss.ext.qlf[subjectHits(ols), ]
  gene.exp <- data.frame(gene.exp$logFC)
  volcano.gp <- unique(i$volcano.gp)
  gene.exp$volcano.gp <- volcano.gp
  gene.exp.dv <- rbind(gene.exp.dv, gene.exp)
}

colnames(gene.exp.dv) <- c("log2FC", "volcano_group")

setwd("~/Thesis/")

## plot results
png(filename='enh2gene_exp.png', height = 700, width=700)
ggplot(gene.exp.dv,
       aes(volcano_group, log2FC)) +
  geom_boxplot(notch = T, outlier.shape = NA, 
               fill = c("red1", "green4", "grey")) +
geom_hline(
    yintercept = 0,
    col = "blue",
    linetype = "dotted",
    size = 1
  ) +
  coord_cartesian(ylim = c(-4, 4)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"))
dev.off()

#### t test
ttest.df <- data.frame()

dv.gps <- list(1,2,4,5,7,8)

for(i in dv.gps){
x <- subset(gene.exp.dv, gene.exp.dv$dv.group == 8)$logFC
y <- subset(gene.exp.dv, gene.exp.dv$dv.group == 9)$logFC

t <- t.test(x, y, alternative = "greater")

t.res <- data.frame(i, t$p.value)
ttest.df <- rbind(ttest.df, t.res)
}
