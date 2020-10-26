library(edgeR)

## for INS vs PNET - switch order of samples
quant_counts <- quant_counts[, c(12:20,1:11)]


## groups 
group.rna <- factor(c(rep("C", 7), rep("T", 11)))


## create DGEList
cds.rna <- DGEList(quant_counts, group = group.rna)

## filter low count reads
cds.rna <- cds.rna[rowSums(1e+06 * cds.rna$counts/expandAsMatrix(cds.rna$samples$lib.size, dim(cds.rna)) > 1) >= 3, ]
# dim( cds )

## re-compute library sizes
cds.rna$samples$lib.size <- colSums(cds.rna$counts)

## calculate normalisation factors
cds.rna <- calcNormFactors(cds.rna)

## Estimating Disperions

# glm
cds.rna <- estimateDisp(cds.rna)


## QLF testing
fit <- glmQLFit(cds.rna)
qlf <- glmQLFTest(fit, coef = 2)
topTags.rna.qlf <- topTags(qlf, n = nrow(cds.rna$counts))
topTags.rna.qlf.df <- topTags.rna.qlf$table


## get counts and boxplot
rna.glm <- data.frame(cpm(cds.rna))





