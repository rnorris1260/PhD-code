### prepare count matrix

setwd("~/R_workspaces/count_mat/")

load("counts.rda")

count_matrix <- counts$counts

colnames(count_matrix) <- gsub(".bam", "", colnames(count_matrix))


group <- factor(c(rep("HI", 6), rep("INS", 12), rep("PNET", 20)))
                             
design <- model.matrix(~group)

library(RUVSeq)
library(edgeR)

####### Two pass normalisation ####### using ruv counts for second pass

### first pass
y.1 <- DGEList(counts=count_matrix, group=group)

y.1 <- calcNormFactors(y.1, method = "upperquartile")

y.1 <- estimateDisp(y.1, design = design)

fit <- glmFit(y.1, design)

res <- residuals(fit, type="deviance")
#We can use all the regions to estimate the factors of unwanted variation.
controls <- rownames(count_matrix)
set.ruv <- RUVr(count_matrix, controls, k=1, res)

enh.RUV <- set.ruv$normalizedCounts


### second pass
## from INS_HI_PNET_nc
y.2 <- DGEList(counts=enh.RUV, group=group)

y.2 <- calcNormFactors(y.2,  method = "RLE")

y.2 <- estimateDisp(y.2)

nc.2 <- cpm(y.2)




