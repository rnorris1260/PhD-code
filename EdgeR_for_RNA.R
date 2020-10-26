library(edgeR)

## for INS vs PNET - switch order of samples
quant_counts <- quant_counts[, c(12:20,1:11)]


## groups 
group.rna <- factor(c(rep("HI", 7), rep("INS", 11)))
group.rna <- factor(c(rep("PNET", 9), rep("INS", 11)))
group.rna <- factor(c(rep("HI", 7), rep("PNET", 9)))


## create DGEList
cds.rna <- DGEList(quant_counts, group = group.rna)

## filter low count reads
# keep <- rowSums(cpm(cds)>2) >= 3
# cds <- cds[keep, , keep.lib.sizes = F]

cds.rna <- cds.rna[rowSums(1e+06 * cds.rna$counts/expandAsMatrix(cds.rna$samples$lib.size, dim(cds.rna)) > 1) >= 3, ]
# dim( cds )

## re-compute library sizes
cds.rna$samples$lib.size <- colSums(cds.rna$counts)

## expression set
#counts.rna <- round(cds.rna$counts,)
#set.rna <- newSeqExpressionSet(as.matrix(counts.rna))

## calculate normalisation factors
cds.rna <- calcNormFactors(cds.rna)

## Estimating Disperions

# glm
#design.rna <- model.matrix(~group.rna)
cds.rna <- estimateDisp(cds.rna)


## QLF testing
fit <- glmQLFit(cds.rna)
qlf <- glmQLFTest(fit, coef = 2)
topTags.rna.qlf <- topTags(qlf, n = nrow(cds.rna$counts))
topTags.rna.qlf.df <- topTags.rna.qlf$table

## exact test
topTags.INS.PNET <- exactTest(cds.rna, pair = c("PNET", "INS"))
topTags.INS.PNET.df <- topTags.INS.PNET$table

exact.rna <- exactTest(cds.rna, pair = c("HI", "INS"))
topTags.rna.INS.HI <- topTags(exact.rna, n = nrow(cds.rna$counts))
topTags.rna.INS.HI.df <- topTags.rna.INS.HI$table

topTags.PNET.HI <- topTags(qlf, n = nrow(cds.rna$counts))
topTags.rna.PNET.HI.df <- topTags.PNET.HI$table


## get counts and boxplot
INS11_HI7.rna.glm <- data.frame(cpm(cds.rna))
boxplot(log10(INS11_HI7.rna.glm), ylim = c(0, 2))



### "Testing"
# Default behavior ("auto") is to use most complex dispersions found in data object.
# de.rna.hg38 <- exactTest(cds.rna, pair = c("HI", "INS"))
# de_RNA.hg38.df <- de.rna.hg38$table
# n <- nrow(de_RNA.hg38.df)
# topTags.rna.hg38 <- topTags(de.rna.hg38, n = n, sort.by = "logFC")
# topTags.rna.hg38 <- topTags.rna.hg38$table





rna.qlf.counts <- cpm(fit)

library("RUVSeq")
### RUVseq
#Again, we can use all the genes to estimate the factors of unwanted variation.

# top <- topTags(qlf, n=nrow(set.rna))$table
# empirical <- rownames(set.rna)[which(!(rownames(set.rna) %in% rownames(top)[1:5000]))]
# # Here, we consider all but the top 5,000 genes as ranked by edgeR p-values.
# set.rna <- RUVg(set.rna, empirical, k=1)
# pData(set.rna)
# plotRLE(set.rna, outline=FALSE, ylim=c(-4, 4), col=colors[x])
# plotPCA(set.rna, col=colors[x], cex=1.2)



### Merge counts with DE results

##prepare data

topTags.rna.INS.HI.df$TX_ID <- rownames(topTags.rna.INS.HI.df)
topTags.rna.qlf.df$TX_ID <- rownames(topTags.rna.qlf.df)


topTags.INS.PNET.df$TX_ID <- rownames(topTags.INS.PNET.df)
INS11_HI7.rna$TX_ID <- rownames(INS11_HI7.rna)

topTags.rna.PNET.HI.df$TX_ID <- rownames(topTags.rna.PNET.HI.df)


                                ### counts      ### all DE res
#rna.res.hg38 <- dplyr::left_join(INS11_HI7.rna, topTags.rna.hg38)


### then merge all res with txdb
#rna.res.hg38.txdb <- dplyr::left_join(rna.res.hg38, tx2gene.38)### select sig 
#rna.res.INS_DE <- subset(rna.res.hg38.txdb, rna.res.hg38.txdb$logFC > 1 & rna.res.hg38.txdb$FDR < 0.05)


#### or join topTags with txdb
rna.res.INS_HI.exact <- dplyr::left_join(topTags.rna.INS.HI.df, tx2gene.38)
rna.res.INS_HI.qlf <- dplyr::left_join(topTags.rna.qlf.df, tx2gene.38)

rna.res.INS_HI <- na.omit(rna.res.INS_HI)
rna.res.INS_HI.qlf <- na.omit(rna.res.INS_HI.qlf)


rna.res.INS_HI <- dplyr::left_join(rna.res.INS_HI, quant_counts.df)

rna.res.INS_PNET <- dplyr::left_join(topTags.INS.PNET.df, tx2gene.38)
rna.res.INS_PNET <- dplyr::left_join(rna.res.INS_PNET, quant_counts.df)

rna.res.INS_PNET <- na.omit(rna.res.INS_PNET)

rna.res.PNET_HI <- dplyr::left_join(topTags.rna.PNET.HI.df, tx2gene.38)

### join with counts
INS11_HI7.rna.glm$TX_ID <- rownames(INS11_HI7.rna.glm)
rna.res.hg38.counts <- dplyr::left_join(rna.res.hg38, INS11_HI7.rna.glm)
rna.res.hg38.counts <- na.omit(rna.res.hg38.counts)


### select only protein-coding transcripts
TX_types <- c("protein_coding", "processed_transcript", "lncRNA", "misc_RNA", "TR_C_gene", "IR_D_gene", "retained_intron",
              "snRNA", "miRNA", "snoRNA", "sRNA", "TR_D_gene", "IG_J_gene", "IG_V_gene", "scRNA")

rna.res.pc.INS_HI <- rna.res.INS_HI[rna.res.INS_HI$TX_biotyoe %in% TX_types, ]
                            



#### generalised linear model

# dge <- DGEList(counts = quant$counts,
#                group = group)
# 
# dge$samples <- merge(dge$samples,
#                      as.data.frame(colData(airway)),
#                      by = 0)
# dge$genes <- data.frame(name = names(rowRanges(airway)),
#                         stringsAsFactors = FALSE)
#dge <- calcNormFactors(dge)



# Next we setup the design matrix and estimate the dispersion (variance). There are multiple ways to do this, and the weird two-step procedure is necessary.
design <- model.matrix(~dge$samples$group)
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

# Now we do a glmFit(), similar to limma
fit <- glmFit(dge, design)

# Now it is time to do a test and extract the top hits
lrt <- glmLRT(fit, coef = 2)
topTags_RNA <- topTags(lrt)


## RNA
topTags_RNA.tab <- topTags_RNA$table
topTags_RNA.tab$logqval <- -log10(topTags_RNA.tab$FDR)

## plot Volcano
setwd("/imppc/labs/lplab/rnorris/Analysis/07_2018_RNA/plots/volcano/")

png(filename='INS10_vs_HI7_RNA_volcano_EdgeR_2.png', height =700, width=700)
par(mar=c(8, 8, 3, 3))
plot(topTags_RNA.tab$logFC, topTags_RNA.tab$logqval, pch=20, 
     main = "", 
     col = "grey", ylab = "", xlab = "",
     cex.axis = 1.5,
     las = 1)  
with(subset(topTags_RNA.tab, FDR<=0.01 & logFC<=-1), points(logFC, logqval, pch=20, col="green"))
with(subset(topTags_RNA.tab, FDR<=0.01 & logFC>=1), points(logFC, logqval, pch=20, col="red"))
title(ylab="RNA -log10 q-value", line=3.5, cex.lab=2)
title(xlab="RNA log2FC (INSvsHI)", line=3.5, cex.lab=2)
dev.off()


## select significant results

topTags_RNA.tab <- topTags_RNA$table
INS_DE.glm.TX <- subset(topTags_RNA.tab, topTags_RNA.tab$logFC >= 1 & topTags_RNA.tab$FDR <= 0.01)
HI_DE.glm.TX <- subset(topTags_RNA.tab, topTags_RNA.tab$logFC <= -1 & topTags_RNA.tab$FDR <= 0.01)




#### get gene names for transcript ids

mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
TX_ids <- as.vector(rownames(topTags_RNA.tab))


result <- getBM(attributes=c("ensembl_transcript_id_version", 'external_gene_name'), 
                filters = "ensembl_transcript_id_version", 
                values = TX_ids, 
                mart= mart)

result$ensembl_transcript_id <- GeneStructureTools::removeVersion(result$ensembl_transcript_id_version)


#### calculate standard deviation for human islet RNA-seq - results with lowest SD more likely to represent true results

rna.res.hg38.txdb$HI_SD <- apply(rna.res.hg38.txdb[, 1:7], 1, sd)

### mean SD of tpm for HI samples is 255.
rna.res.SDlim <- subset(rna.res.hg38.txdb, rna.res.hg38.txdb$HI_SD <= 255)

### now lets get all the INS.DE genes for the most consistent HI RNA-seq results
rna.res.SDlim.INS <- subset(rna.res.SDlim, rna.res.SDlim$logFC > 1 & rna.res.SDlim$FDR < 0.05)
genes.SDlim.INS <- unique(data.frame(rna.res.SDlim.INS$gene_name))


### lets assume for a minute that the most interesting results will come from overexpression of genes that are already expressed in islets - not de novo expression
### so lets subset the rna results for genes that have a minimal level of expression in islets

rna.res.SDlim.INS$HI_mean <- apply(rna.res.SDlim.INS[, 1:7], 1, mean)
rna.res.SDlim.HIex <- subset(rna.res.SDlim.INS, rna.res.SDlim.INS$HI_mean >= 5)
genes.SDlim.HIex <- unique(data.frame(rna.res.SDlim.HIex$gene_name))
write.csv(genes.SDlim.HIex, file = "genes.SDlim.HIex.csv")



###### Gene Set Enrich,ent Analysis with SetRank

library("SetRank")

annotation.tab <- tx2gene.38[, c(2,3,5)]

buildSetCollection(annotation.tab)
