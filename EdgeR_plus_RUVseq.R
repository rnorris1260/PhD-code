### prepare count matrix

setwd("~/R_workspaces/count_mat/")

save(count_matrix, file = "INS_HI_PNET_INSlike_counts.rda")

load("INS_HI_PNET_beta_counts.rda")

count_matrix <- counts$counts

colnames(count_matrix) <- gsub(".bam", "", colnames(count_matrix))

## select HI and non-beta PNETs
count_matrix_non_beta <- count_matrix[, c(1:6,19:24,26:28,30,32:37)]
count_matrix_non_beta <- count_matrix_non_beta[, -c(15:21)]

#count_matrix.2b <- count_matrix.2[, c(1:18,29,31,25,38)]

group <- factor(c(rep("HI", 6), rep("INS", 12), rep("PNET", 20)))

group <- factor(c(rep("HI", 6), rep("INS", 12), rep("PNET", 3)))

group.nb <- factor(c(rep("HI", 6), rep("PNET", 9)))

#set <- newSeqExpressionSet(as.matrix(count_matrix.2))
                             
design <- model.matrix(~group)


#design <- model.matrix(~group.2, data=pData(set))


library(RUVSeq)
library(edgeR)

####### Two pass normalisation ####### using ruv counts for second pass

### first pass
y.1 <- DGEList(counts=count_matrix, group=group)

y.1 <- calcNormFactors(y.1, method = "upperquartile")

y.1 <- estimateDisp(y.1, design = design)

fit <- glmFit(y.1, design)

res <- residuals(fit, type="deviance")
#Again, we can use all the genes to estimate the factors of unwanted variation.
controls <- rownames(count_matrix)
set.ruv <- RUVr(count_matrix, controls, k=1, res)

INS_HI_PNET.RUV <- set.ruv$normalizedCounts


### second pass
## from INS_HI_PNET_nc
y.2 <- DGEList(counts=INS_HI_PNET.RUV, group=group)

y.2 <- calcNormFactors(y.2,  method = "RLE")

y.2 <- estimateDisp(y.2)

nc.2 <- cpm(y.2)




## DE

# INS vs HI
INSvsHI <- exactTest(y.2, pair = c("HI", "INS"))
topTags.INS.HI.enh <- topTags(INSvsHI, n=nrow(count_matrix))
topTags.INS.HI.enh.df <- topTags.INS.HI.enh$table

# INS vs PNET
INSvsPNET <- exactTest(y.2, pair = c("PNET", "INS"))
topTags.INS.PNET <- topTags(INSvsPNET, n=nrow(count_matrix))
topTags.INS.PNET.df <- topTags.INS.PNET$table

# PNET vs HI
PNETvsHI <- exactTest(y.2, pair = c("HI", "PNET"))
topTags.PNET.HI <- topTags(PNETvsHI, n=nrow(count_matrix))
topTags.PNET.HI.df <- topTags.PNET.HI$table

topTags.PNET.HI.df[, 5:7] <- region.list(topTags.PNET.HI.df)
topTags.PNET.HI.df <- topTags.PNET.HI.df[, c(5:7,1,3,4)]
rownames(topTags.PNET.HI.df) <- c(1:nrow(topTags.PNET.HI.df))
PNET_HI_DA <- subset(topTags.PNET.HI.df, topTags.PNET.HI.df$logFC > 1.5 &
                             topTags.PNET.HI.df$FDR < 0.05)


### merge data frames
INS.HI.dv <- topTags.INS.HI.enh.df[, c(1,4)]
INS.HI.dv[, 3:5] <- region.list(INS.HI.dv)
INS.HI.dv <- INS.HI.dv[, c(3:5,1,2)]
rownames(INS.HI.dv) <- c(1:nrow(INS.HI.dv))
colnames(INS.HI.dv)[c(4:5)] <- c("HI.logFC", "HI.FDR")

INS.PNET.dv <- topTags.INS.PNET.df[, c(1,4)]
INS.PNET.dv[, 3:5] <- region.list(INS.PNET.dv)
INS.PNET.dv <- INS.PNET.dv[, c(3:5,1,2)]
rownames(INS.PNET.dv) <- c(1:nrow(INS.PNET.dv))
colnames(INS.PNET.dv)[c(4:5)] <- c("PNET.logFC", "PNET.FDR")


dv.df <- dplyr::left_join(INS.HI.dv, INS.PNET.dv)



INS_HI_PNET.cor <- cor(nc.2, method = "spearman")



library(ComplexHeatmap)

png(file='INS_HI_PNET_corS.plusRUV_2pass.png', height = 700, width = 900)

Heatmap(INS_HI_PNET.cor, clustering_method_columns = method,
        clustering_method_rows = method, col = col, 
        right_annotation=ra, bottom_annotation = ca,
        heatmap_legend_param = list(title="Spearman correlation coefficient"),
        column_dend_height = unit(3, "cm"),
        row_dend_width = unit(3, "cm"),
        row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10))

dev.off()


png(filename='INS_HI_PNET_post-norm.png', height = 700, width=1000)
par(mar = c(10, 4, 4, 2))

boxplot(log2(nc.2 +1), col = "powderblue", 
        main = "HI (x6), INS (x12), PNET (x20) counts post-norm", 
        ylab = "log2 counts post-norm",
        las = 2)

dev.off()



### PCA
#DESeq2
library(DESeq2)
counts.norm <- as.matrix(round(nc.2, 0))

col_data <- data.frame("tissue"=colnames(counts.norm))
col_data$cell_type <- c(rep("HI", 6), rep("INS", 12), rep("PNET", 20))


col_data[c(19,20,23,24,26,27,21,22,37), 2] <- "PNET_ARX"
col_data[c(29,31,25,38), 2] <- "PNET_PDX"
col_data[c(33,35,34,28,36,30,32), 2] <- "PNET_DP"


col_data$cell_type <- factor(x = col_data$cell_type, levels = c("HI", "INS", "PNET_ARX", "PNET_PDX", "PNET_DP"))


dds <- DESeqDataSetFromMatrix(countData=counts.norm, 
                              colData=col_data,
                              design = ~cell_type)


rld <- rlogTransformation(dds, blind=TRUE)
rld.vsd <- vsd(dds, blind=TRUE)

plotPCA(dds, intgroup = "cell_type", ntop=2000) +
        geom_text(aes(label=name),vjust=2,check_overlap = TRUE, size = 4)





### direct from tutorial

z <- DGEList(count_matrix.2, group=group.2)
z <- calcNormFactors(y, method="upperquartile") #RUVSeq: Remove Unwanted Variation from RNA-Seq Data

z <- estimateDisp(y, design = design.2)

fit.z <-glmFit(z, design.2)

res.z <- residuals(fit.z, type="deviance")
#Again, we can use all the genes to estimate the factors of unwanted variation.
set.z <- RUVr(count_matrix.2, genes, k=1, res)
z.cn <- set.z$normalizedCounts
z.cn.log <- log2(z.cn +1)



###### Do enhancers correspond to changes in gene expression on a general level?

###### dv group 8 vs all rna res

rownames(tss.ext.qlf) <- c(1:nrow(tss.ext.qlf))
tss.ext.qlf.gr <- toGRanges(tss.ext.qlf[, c(9,14,15,1)])

rna.vs.df8 <- findOverlaps(dv.gp8, tss.ext.qlf.gr)

rna.res.gp8 <- data.frame(tss.ext.qlf[subjectHits(rna.vs.df8), ])
rna.res.gp8 <- rna.res.gp8[, c(1,8)]




### non-beta

####### Two pass normalisation ####### using ruv counts for second pass

### first pass
y.1.nb <- DGEList(counts=count_matrix_non_beta, group=group.nb)

y.1.nb <- calcNormFactors(y.1.nb, method = "upperquartile")

design.nb <- model.matrix(~group.nb)

y.1.nb <- estimateDisp(y.1.nb, design = design.nb)

fit.nb <- glmFit(y.1.nb, design.nb)

res.nb <- residuals(fit.nb, type="deviance")
#Again, we can use all the genes to estimate the factors of unwanted variation.
controls <- rownames(count_matrix_non_beta)
set.ruv.nb <- RUVr(count_matrix_non_beta, controls, k=1, res.nb)

HI_PNET.RUV <- set.ruv.nb$normalizedCounts


### second pass
## from INS_HI_PNET_nc
y.2.nb <- DGEList(counts=HI_PNET.RUV, group=group.nb)

y.2.nb <- calcNormFactors(y.2.nb,  method = "RLE")

y.2.nb <- estimateDisp(y.2.nb)

nc.nb <- cpm(y.2.nb)

# fit <- glmFit(y.2)



## DE

# INS vs HI
PNETvsHI <- exactTest(y.2.nb, pair = c("HI", "PNET"))
topTags.PNETnb.HI <- topTags(PNETvsHI, n=nrow(count_matrix_non_beta))
topTags.PNETnb.HI.df <- topTags.PNETnb.HI$table

topTags.PNETnb.HI.df[, 5:7] <- region.list(topTags.PNETnb.HI.df)
topTags.PNETnb.HI.df <- topTags.PNETnb.HI.df[, c(5:7,1:4)]

PNETnb.DA <- subset(topTags.PNETnb.HI.df, topTags.PNETnb.HI.df$logFC > 1.5 &
                            topTags.PNETnb.HI.df$FDR < 0.05)



# y.1 <- DGEList(counts=counts(set), group=group.2)
# y.1 <- calcNormFactors(y.1, method="upperquartile") #RUVSeq: Remove Unwanted Variation from RNA-Seq Data
# 
# y.1 <- estimateDisp(y.1, design.2)
# 
# fit <- glmFit(y.1, design.2)
# 
# # RUV
# res <- residuals(fit, type="deviance")
# #Again, we can use all the genes to estimate the factors of unwanted variation.
# #set <- betweenLaneNormalization(set, which="upper")
# controls <- rownames(count_matrix.2)
# set <- RUVr(set, controls, k=1, res)
# pData(set.1)
# 
# #W_1 <- set.1$W_1
# 
# #design..2 <- model.matrix(~group.2 + W_1, data=pData(set.1))
# y.2 <- DGEList(counts=counts(set), group=group.2)
# y.2 <- calcNormFactors(y, method = "upperquartile")
# y.2 <- estimateDisp(y, design..2)
# 
# nc.1 <- cpm(y)
# 
# fit <- glmFit(y, design)
# lrt <- glmLRT(fit, coef=2)
# topTags(lrt)
