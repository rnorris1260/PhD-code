##### fast gene set enrichment analysis with FGSEA

library(fgsea)
library(ggplot2)

## get largest logFC for each gene
gene.list.mod <- function(gene_list){
  gene_list <- gene_list[order(gene_list$logFC, decreasing = T),]
  gene_list <- gene_list[!duplicated(gene_list$gene_name), ]
  gene_list <- na.omit(gene_list)
}


## RNA blacklist
TX_biotype <- data.frame(table(rna.res.INS_HI$TX_biotype)[rna.res.INS_HI$TX_biotype])
TX_biotype <- TX_biotype[!duplicated(TX_biotype$Var1), ]

blacklist <- c("nonsense_mediated_decay", "non_stop_decay", "TEC")

#rna.res.INS_HI <- na.omit(rna.res.INS_HI)
rna.res.INS_HI.qlf <- na.omit(rna.res.INS_HI.qlf)
rna.res.INS_HI.exact <- na.omit(rna.res.INS_HI.exact)

#rna.res.INS_HI <- rna.res.INS_HI[rna.res.INS_HI$TX_biotyoe %ni% blacklist, ]
rna.res.INS_HI.qlf <- rna.res.INS_HI.qlf[rna.res.INS_HI.qlf$TX_biotype %ni% blacklist, ]
rna.res.INS_HI.exact <- rna.res.INS_HI.exact[rna.res.INS_HI.exact$TX_biotype %ni% blacklist, ]

#all.genes <- rna.res.INS_HI[, c(6,1)]
all.genes.qlf <- rna.res.INS_HI.qlf.canon[, c(8,1)]
all.genes.exact <- rna.res.INS_HI.exact[, c(7,1)]

all.genes <- gene.list.mod(all.genes.qlf)

## create ranks 
ranks <- tibble::deframe(all.genes)



#### GSEA code
GSEA = function(gene_list, GO_file, pval) {
  set.seed(54321)
  library(dplyr)
  library(gage)
  library(fgsea)
  
  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
  myGO = fgsea::gmtPathways(GO_file)
  
  fgRes <- fgsea::fgseaMultilevel(pathways = myGO, 
                        stats = gene_list,
                        minSize=15,
                        maxSize=500,
                        ) %>% 
    as.data.frame() %>% 
    dplyr::filter(padj < !!pval)
  #print(dim(fgRes))
  
  ## Filter FGSEA by using gage results. Must be significant and in same direction to keep 
  gaRes = gage::gage(gene_list, gsets=myGO, same.dir=TRUE, set.size =c(15,500))
  
  ups = as.data.frame(gaRes$greater) %>% 
    tibble::rownames_to_column("Pathway") %>% 
    dplyr::filter(!is.na(p.geomean) & q.val < pval ) %>%
    dplyr::select("Pathway")
  
  downs = as.data.frame(gaRes$less) %>% 
    tibble::rownames_to_column("Pathway") %>% 
    dplyr::filter(!is.na(p.geomean) & q.val < pval ) %>%
    dplyr::select("Pathway")
  
  #print(dim(rbind(ups,downs)))
  keepups = fgRes[fgRes$NES > 0 & !is.na(match(fgRes$pathway, ups$Pathway)), ]
  keepdowns = fgRes[fgRes$NES < 0 & !is.na(match(fgRes$pathway, downs$Pathway)), ]
  
  ### sort & subset res
  keepups <- keepups[order(keepups$pval), ]
  keepups <- keepups[keepups$padj < 0.05, ]
  
  ### Collapse redundant pathways
  Up <- collapsePathways(keepups,
                           pathways = myGO, stats = gene_list,
                           nperm = 1000)

  #Up = fgsea::collapsePathways(keepups, pathways = myGO, stats = gene_list,  nperm = 1000)
  Down = fgsea::collapsePathways(keepdowns, myGO, gene_list,  nperm = 1000)
  
  fgRes = fgRes[ !is.na(match(fgRes$pathway, 
                              c( Up$mainPathways, Down$mainPathways))), ] %>% 
    arrange(desc(NES))
  fgRes$pathway = stringr::str_replace(fgRes$pathway, "REACTOME_" , "")
  
  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  #filtRes = rbind(head(fgRes, n = 10),
      #            tail(fgRes, n = 10 ))
  return(fgRes)
}
  


##### do GSEA
#gsea.res.exact <- GSEA(gene_list = ranks, GO_file = "~/pathways/gmt_files/c2.cp.reactome.v7.1.symbols.gmt", pval = 0.05)

gsea.res.qlf <- GSEA(gene_list = ranks, GO_file = "~/pathways/gmt_files/c2.cp.reactome.v7.1.symbols.gmt", pval = 0.05)
gsea.res.qlf.check <- GSEA(gene_list = ranks, GO_file = "~/pathways/gmt_files/c2.cp.reactome.v7.1.symbols.gmt", pval = 0.05)
gsea.res.qlf.canon <- GSEA(gene_list = ranks, GO_file = "~/pathways/gmt_files/c2.cp.reactome.v7.1.symbols.gmt", pval = 0.05)
gsea.res.qlf.canon.500 <- GSEA(gene_list = ranks, GO_file = "~/pathways/gmt_files/c2.cp.reactome.v7.1.symbols.gmt", pval = 0.05)


g = ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
  geom_segment( aes(reorder(pathway, NES), xend=pathway, y=0, yend=NES)) +
  geom_point( size=5, aes( fill = Enrichment),
              shape=21, stroke=2) +
  scale_fill_manual(values = c("Down-regulated" = "dodgerblue",
                               "Up-regulated" = "firebrick") ) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GSEA - Biological Processes") + 
  theme_minimal()

output = list("Results" = fgRes, "Plot" = g)
return(output)


#NEW GSEA code

### select pathways associated with insulinoma-specific enhancers
## plus chromatin modifiers

gsea.res.select <- gsea.res.qlf[c(2,3,9:11,13,15,16,18,19), ]

# FC 1.5
GSEA.vs.DI.all <- data.frame(gsea.res.select[, 1])
GSEA.vs.DI.b <- data.frame(gsea.res.select[, 1])
GSEA.vs.DI.all.FC1 <- data.frame(gsea.res.select[, 1])
GSEA.vs.DI.all.2 <- data.frame(gsea.res.qlf.check[, 1])

GSEA.vs.DI.all.canon <- data.frame(gsea.res.qlf.canon.500[, 1])
GSEA.vs.DI.all.canon.chk <- data.frame(gsea.res.select[, 1])
GSEA.vs.DI.all.canon.b <- data.frame(gsea.res.select[, 1])

### select pathways to test for enhancer association
colnames(GSEA.vs.DI.all.canon)[1] <- c("pathway")
colnames(GSEA.vs.DI.all.canon.chk)[1] <- c("pathway")
colnames(GSEA.vs.DI.all.canon.b)[1] <- c("pathway")

## select genes from gsea pathways
function(gsea, r){
  genes <- data.frame(gsea[r,8]) 
  colnames(genes) <- "genes"
  rna.res <- tss.ext.qlf.cut[tss.ext.qlf.cut$gene_name %in% genes$genes, ]
  rna.res <- rna.res[order(rna.res$logFC, decreasing = T), ]
  rna.res.genes <- rna.res[!duplicated(rna.res$gene_name), ]
}

library(regioneR)

genome <- getGenome("hg38")

for(i in 1:9){
  genes <- gsea.genes(gsea.res.select, i)
  genes.gr <- regioneR::toGRanges(genes)
  
  pt <- overlapPermTest(A=genes.gr, B=toGRanges(INS.specSE.enh), ntimes=1000, genome = genome,
                        alternative = "greater", verbose = TRUE)
  GSEA.vs.DI.all.canon.chk[i,18] <- pt$numOverlaps$pval
}

colnames(GSEA.vs.DI.all.canon.chk)[18] <- c("INS.SE.pval")

GSEA.vs.DI.all.canon.chk$INS.SE.FDR <- p.adjust(GSEA.vs.DI.all.canon.chk$INS.SE.pval, method = "bonferroni", n = 9)


GSEA.vs.DI.all.FDR <- GSEA.vs.DI.all.canon[, c(1,3,5,7,9,11,13,15,17,19)]
GSEA.vs.DI.all.FDR <- GSEA.vs.DI.all.canon.chk[, c(1,3,5,7,9,11,13,15,17,19)]
GSEA.vs.DI.b.FDR <- GSEA.vs.DI.b.canon[, c(1,3,5,7,9,11,13,15,17,19)]


### remove pathways that have no overlap with enhancers
gsea.res.select <- gsea.res.qlf.canon.500[c(2:4,9:11,14:17), ]
GSEA.vs.DI.all.canon <- GSEA.vs.DI.all.canon[c(2:4,9:11,14:17), ]

gsea_investig <- dplyr::left_join(gsea.res.qlf, GSEA.vs.DI.all)

### results of overlap of gsea pathways with super-enhancers suggest there are 
### chromatin modifying enzymes outside of the HAT and HDM groups with novel super-enhancer(s)


### HATs
HAT.genes <- data.frame(gsea.res.select[1,8])
colnames(HAT.genes) <- "gene_name"


HDM.genes <- data.frame(gsea.res.select[2,8])
colnames(HDM.genes) <- "gene_name"


chrom.mod.genes <- data.frame(gsea.res.select[3,8])
colnames(chrom.mod.genes) <- "gene_name"

chrom.mod.extra <- subset(chrom.mod.genes, chrom.mod.genes$gene_name %ni% HAT.genes$gene_name &
                                            chrom.mod.genes$gene_name %ni% HDM.genes$gene_name)

chrom.mod.extra.res <- rna.res.INS_HI.qlf[rna.res.INS_HI.qlf$gene_name %in% chrom.mod.extra$gene_name, ]

chrom.mod.tss.ext <- tss.ext.qlf[tss.ext.qlf$gene_name %in% chrom.mod.extra.res$gene_name, ]

rownames(chrom.mod.tss.ext) <- c(1:nrow(chrom.mod.tss.ext))

INS.SE.inCM <- findOverlaps(toGRanges(INS.spec.SEs.enh), toGRanges(chrom.mod.tss.ext[, c(9,14,15)]))

INS.SE.CM.genes <- chrom.mod.tss.ext[rownames(chrom.mod.tss.ext) %in% subjectHits(INS.SE.inCM), ]

INS.SE.CM.genes <- INS.SE.CM.genes[!duplicated(INS.SE.CM.genes$TX_ID), ]


## merge INS.SE with chrom mod extra to see how many SE hits there are for each gene

INS.spec.SEs.enh <- toGRanges(INS.spec.SEs.enh) 

chrom.mod.tss.ext <- toGRanges(chrom.mod.tss.ext[, c(9,14,15, 1:8,10:13)])

INS.SE.CM.res <- cbind(data.frame(INS.spec.SEs.enh[queryHits(INS.SE.inCM),]),
                                       data.frame(chrom.mod.tss.ext[subjectHits(INS.SE.inCM),]))

INS.SE.CM.vsgp1 <- findOverlaps(toGRanges(INS.SE.CM.res), toGRanges(df1b))
INS.SE.CM.vsgp1 <- INS.SE.CM.res[rownames(INS.SE.CM.res) %in% queryHits(INS.SE.CM.vsgp1), ]

INS.SE.CM.vsgp8 <- findOverlaps(toGRanges(INS.SE.CM.res), toGRanges(df8b))
INS.SE.CM.vsgp8 <- INS.SE.CM.res[rownames(INS.SE.CM.res) %in% queryHits(INS.SE.CM.vsgp8), ]



# select only up-regulated transcripts
INS.SE.CM.res.upreg <- subset(INS.SE.CM.res, INS.SE.CM.res$logFC > 1)

INS.SE.CM.res.gene.rep <- data.frame(table(INS.SE.CM.res.upreg$gene_name)[INS.SE.CM.res.upreg$gene_name])
INS.SE.CM.res.gene.rep <- INS.SE.CM.res.gene.rep[!duplicated(INS.SE.CM.res.gene.rep$Var1), ]

### load INS-specfic super-enhancers
setwd("~/R_workspaces/objects_160620/")
load("INS.spec.SEs.enh.rda")
INS.spec.SEs.enh <- data.frame(INS.spec.SEs.enh)



### compare genes from pathways associated with different enhancer groups
### the two pathways the same results in vsall and vsbeta are nervous system dev and RTK signalling

tss.ext.qlf.cut <- tss.ext.qlf[, c(9,16,17,1,5:8,10:15)]

## regular enhancer data 
INS_REs <- read.bed("INS_REs.bed")
colnames(INS_REs)[c(4,5)] <- c("enh_id", "SI")

## merge regular enhancers with SEs
rownames(df8.all) <- c(1:nrow(df8.all))
df8.enh.ols <- findOverlaps(toGRanges(INS_REs), toGRanges(df8.all))

df8.enh <- cbind(data.frame(INS_REs[queryHits(df8.enh.ols), ]),
                        data.frame(df8.all[subjectHits(df8.enh.ols), ]))
rownames(df8.enh) <- c(1:nrow(df8.enh))


enh.ol.gsea <- function(gs, enh_gp){
  rownames(enh_gp) <- c(1:nrow(enh_gp))
  genes <- data.frame(gsea.res.select[gs,8])
  colnames(genes) <- "gene_name"
  rna <- tss.ext.qlf.cut[tss.ext.qlf.cut$gene_name %in% genes$gene_name, ]
  rna <- rna[order(rna$logFC, decreasing = T), ]
  rna <- rna[!duplicated(rna$gene_name), ]
  rownames(rna) <- c(1:nrow(rna))
  ols <- findOverlaps(toGRanges(enh_gp), toGRanges(rna))
  rna.res <- data.frame(rna[rownames(rna) %in% subjectHits(ols), ])
  enh.res <- data.frame(enh_gp[rownames(enh_gp) %in% queryHits(ols), ])
  all.res <- cbind(data.frame(enh_gp[queryHits(ols), ]),
                   data.frame(rna[subjectHits(ols), ]))
  return(list("rna" = rna.res, "enh" = enh.res, "all" = all.res))
}

# pathways:
# 1       RMTS_METHYLATE_HISTONE_ARGININES
# 2                HATS_ACETYLATE_HISTONES
# 3              HDMS_DEMETHYLATE_HISTONES
# 4                      SIGNALING_BY_VEGF
# 5            CHROMATIN_MODIFYING_ENZYMES
# 6                        NEURONAL_SYSTEM
# 7               SIGNALING_BY_RHO_GTPASES
# 8             NERVOUS_SYSTEM_DEVELOPMENT
# 9 SIGNALING_BY_RECEPTOR_TYROSINE_KINASES


get.pathway <- function(p){
  genes <- data.frame(gsea.res.select[p,8])
  colnames(genes) <- "gene_name"
  genes.pathway <- rna.res.INS_HI.qlf.canon[rna.res.INS_HI.qlf.canon$gene_name %in% genes$gene_name, ]
  # genes.logFC <- genes.pathway[, c(8,1,5)]
  # genes.logFC <- genes.logFC[order(genes.logFC$gene_name), ]
  # genes.logFC$gene_name <- factor(genes.logFC$gene_name) 
  return(genes.pathway)
}


## NS dev
NS.dev.all <- enh.ol.gsea(8, df8.enh)
NS.dev.b <- enh.ol.gsea(8, df8b)
NS.dev.all.genes <- NS.dev.all$rna
NS.dev.all2 <- NS.dev.all$all

NS.dev.all2.pc <- subset(NS.dev.all2, NS.dev.all2$TX_biotype == "protein_coding")
NS.dev.all2.pc <- subset(NS.dev.all2.pc, NS.dev.all2.pc$type == "DI-2")
NS.dev.all2.pc <- subset(NS.dev.all2.pc, NS.dev.all2.pc$FDR < 0.05)
NS.dev.all2.pc$gene.freq <- table(NS.dev.all2.pc$gene_name)[NS.dev.all2.pc$gene_name]

NS.dev.LE <- data.frame(gsea.res.select[8,8])

## RTK
RTK.all <- enh.ol.gsea(9, df8.enh)
RTK.b <- enh.ol.gsea(9, df8b)
RTK.all2 <- RTK.all$all

nrow(RTK.LE[!duplicated(RTK.LE$gene_name), ])

RTK.all2.pc <- subset(RTK.all2, RTK.all2$TX_biotype == "protein_coding")
RTK.all2.pc <- subset(RTK.all2.pc, RTK.all2.pc$type == "DI-2")
RTK.all2.pc <- subset(RTK.all2.pc, RTK.all2.pc$FDR < 0.05)
RTK.all2.pc$gene.freq <- table(RTK.all2.pc$gene_name)[RTK.all2.pc$gene_name]


RTK.LE <- data.frame(gsea.res.select[9,8])

## RHO GTPases
RHO_GTP.LE <- get.pathway(7)

nrow(RHO_GTP.LE[!duplicated(RHO_GTP.LE$gene_name), ])

RHO_GTP.all <- enh.ol.gsea(7, df8.enh)
RHO_GTP.b <- enh.ol.gsea(7, df8b)
RHO_GTP.all2 <- RHO_GTP.all$all

RHO_GTP.all2.pc <- subset(RHO_GTP.all2, RHO_GTP.all2$TX_biotype == "protein_coding")
RTK.all2.pc <- subset(RTK.all2.pc, RTK.all2.pc$type == "DI-2")
RHO_GTP.all2.pc <- subset(RHO_GTP.all2.pc, RHO_GTP.all2.pc$FDR < 0.05)
RHO_GTP.all2.pc$gene.freq <- table(RHO_GTP.all2.pc$gene_name)[RHO_GTP.all2.pc$gene_name]

# TFs
RHO_GTP.TFs <- RHO_GTP.all2.pc[RHO_GTP.all2.pc$gene_name %in% TF.list$TF, ]

## SE analysis
RHO_GTP.SE2.SI <- subset(RHO_GTP.SE2.pc, RHO_GTP.SE2.pc$SI >= 3)
RHO_SE_gene_freq <- data.frame(table(RHO_GTP.SE2.SI$gene_name)[RHO_GTP.SE2.SI$gene_name])
RHO_SE_gene_freq <- RHO_SE_gene_freq[!duplicated(RHO_SE_gene_freq$Var1), ]

## VEGF
VEGF <- get.pathway(4)

nrow(VEGF[!duplicated(VEGF$gene_name), ])

VEGF.all <- enh.ol.gsea(9, df8.enh)
VEGF.b <- enh.ol.gsea(9, df8b)
VEGF.all2 <- VEGF.all$all

VEGF.all2.pc <- subset(VEGF.all2, VEGF.all2$TX_biotype == "protein_coding")
VEGF.all2.pc <- subset(VEGF.all2, VEGF.all2$type == "DI-2")
VEGF.all2.pc <- subset(VEGF.all2.pc, VEGF.all2.pc$FDR < 0.05)
VEGF.all2.pc$gene.freq <- table(VEGF.all2.pc$gene_name)[VEGF.all2.pc$gene_name]


## find common genes
common.genes <- NS.dev.all.genes[NS.dev.all.genes$TX_ID %in% RTK.all.genes$TX_ID, ]
common.genes <- common.genes[common.genes$TX_ID %in% RHO_GTP.all2$TX_ID, ]
common.genes <- common.genes[common.genes$TX_ID %in% VEGF.all2$TX_ID, ]

## chrom mod
chrom.mod.all <- enh.ol.gsea(5, df8.enh)
chrom.mod.b <- enh.ol.gsea(5, df8b)
chrom.mod.all.genes <- chrom.mod.all$rna
chrom.mod.all.genes <- subset(chrom.mod.all.genes, chrom.mod.all.genes$FDR < 0.05)
chrom.mod.all2 <- chrom.mod.all$all

chrom.mod.all2.pc <- subset(chrom.mod.all2, chrom.mod.all2$TX_biotype == "protein_coding")
chrom.mod.all2.pc <- subset(chrom.mod.all2.pc, chrom.mod.all2.pc$type == "DI-2")
chrom.mod.all2.pc <- subset(chrom.mod.all2.pc, chrom.mod.all2.pc$FDR < 0.05)
chrom.mod.all2.pc$gene.freq <- table(chrom.mod.all2.pc$gene_name)[chrom.mod.all2.pc$gene_name]

chrom.mod.TFs <- chrom.mod.all2.pc[chrom.mod.all2.pc$gene_name %in% TF.list$TF, ]

chrom.mod.LE <- get.pathway(5)
## chrom mod genes with no associated enhancer activation
chrom.mod.LE.mut <- chrom.mod.LE[chrom.mod.LE$TX_ID %ni% chrom.mod.all2$TX_ID, ]
chrom.mod.LE.mut <- subset(chrom.mod.LE.mut, chrom.mod.LE.mut$TX_biotype == "protein_coding")
chrom.mod.LE.mut <- subset(chrom.mod.LE.mut, chrom.mod.LE.mut$logFC > 1 
                           & chrom.mod.LE.mut$FDR < 0.05)

chrom.mod.noenh <- data.frame(chrom.mod.LE[chrom.mod.LE$gene_name %ni% chrom.mod.all.genes$gene_name, ])
colnames(chrom.mod.noenh) <- "gene_name"
chrom.mod.noenh <- rna.res.INS_HI.qlf.canon[rna.res.INS_HI.qlf.canon$gene_name %in% chrom.mod.noenh$gene_name, ]
chrom.mod.noenh <- subset(chrom.mod.noenh, chrom.mod.noenh$FDR < 0.05 &
                            chrom.mod.noenh$TX_biotype == "protein_coding")

chrom.mod.noenh.SE <- chrom.mod.SE[chrom.mod.SE$gene_name %in% chrom.mod.noenh$gene_name, ]


#### CREs with mutations
mutations <- read.csv(file = "mutations.csv", header = T)
mutations$chr <- paste0("chr", mutations$chr)
rownames(chrom.mod.LE.mut) <- 1:nrow(chrom.mod.LE.mut)
mut.dist <- distanceToNearest(toGRanges(mutations), toGRanges(rna.res.INS_HI.qlf.canon[, c(9,13,14)]))
mut.dist.df <- data.frame(mut.dist)

mut.dist.df.TSS <- subset(mut.dist.df, mut.dist.df$distance == 0)

mut.TSS <- rna.res.INS_HI.qlf.canon[rownames(rna.res.INS_HI.qlf.canon) %in% mut.dist.df.TSS$subjectHits, ]
mut.all <- rna.res.INS_HI.qlf.canon[rownames(rna.res.INS_HI.qlf.canon) %in% mut.dist.df$subjectHits, ]
mut.all <- subset(mut.all, mut.all$TX_biotype == "protein_coding")

all.enh <- topTags.INS.HI.enh.df
all.enh[, 5:7] <- region.list(all.enh)
all.enh <- all.enh[, c(5:7,1:4)]
rownames(all.enh) <- c(1:nrow(all.enh))

enh.ol.mut <- function(genes, enh_gp){
  ols <- findOverlaps(toGRanges(enh_gp), toGRanges(genes[, c(9,13,14)]))
  rna.res <- data.frame(genes[rownames(genes) %in% subjectHits(ols), ])
  enh.res <- data.frame(enh_gp[rownames(enh_gp) %in% queryHits(ols), ])
  all.res <- cbind(data.frame(enh_gp[queryHits(ols), ]),
                   data.frame(genes[subjectHits(ols), ]))
  return(list("rna" = rna.res, "enh" = enh.res, "all" = all.res))
}

CRE_mut <- enh.ol.mut(mut.all, all.enh)

## subset for sharing index
CRE_mut <- subset(CRE_mut, CRE_mut$SI >= 3)



chrom.mod.INSspecSE <- enh.ol.gsea(5, INS.specSE.enh)
HAT.INSspecSE <- enh.ol.gsea(2, INS.specSE.enh)

chrom.mod.vs.SE <- findOverlaps(toGRanges(chrom.mod.INSspecSE$rna), toGRanges(INS.specSE.enh))

chrom.mod.SE <- cbind(data.frame(chrom.mod.INSspecSE$rna[queryHits(chrom.mod.vs.SE), ]),
                      data.frame(INS.specSE.enh[subjectHits(chrom.mod.vs.SE), ]))

chrom.mod.SE.sig <- subset(chrom.mod.SE, chrom.mod.SE$FDR < 0.05)
chrom.mod.SE.sig.cut <- chrom.mod.SE.sig[, c(15:19,4:11)]
chrom.mod.SE.sig.cut$Freq <- table(chrom.mod.SE.sig.cut$gene_name)[chrom.mod.SE.sig.cut$gene_name]

chrom.mod.SE.sig.cut.enh <- chrom.mod.SE.sig.cut[!duplicated(chrom.mod.SE.sig.cut$enh_id), ]

HAT <- enh.ol.gsea(2, df8.all)
HAT.genes <- HAT$rna
HAT.genes <- subset(HAT.genes, HAT.genes$FDR < 0.05 & HAT.genes$TX_biotype == "protein_coding")

### select only protein-coding TX, sharing index > 3, FDR < 0.05
chrom.mod.SE.sig.cut <- subset(chrom.mod.SE.sig.cut, chrom.mod.SE.sig.cut$SI >=3 & 
                               chrom.mod.SE.sig.cut$FDR < 0.05 &
                                 chrom.mod.SE.sig.cut$TX_biotype == "protein_coding")

RMTS <- enh.ol.gsea(1, df8.all)
RMTS.genes <- RMTS$rna


## Neuronal system
neuronal.sys <- enh.ol.gsea(6, df8.all)
neuronal.sys.genes <- neuronal.sys$rna


## enhancers
chrom.mod.all.enh <- chrom.mod.all$enh


NS.dev.all.enh <- NS.dev.all$enh


### MOZ/MORF
MOZ_MORF <- list("BRD1", "BRPF1", "BRPF3", "KAT6A", "KAT6B", "MEAF6")
MOZ_MORF <- rna.res.INS_HI.qlf.canon[rna.res.INS_HI.qlf.canon$gene_name %in% MOZ_MORF, ]
MOZ_MORF <- subset(MOZ_MORF, MOZ_MORF$FDR < 0.05 &
                     MOZ_MORF$TX_biotype == "protein_coding")

### get ratio of enhancers to genes for each gsea pathway
get.ratio <- function(path.num, enh_gp){
  pathway.res <- enh.ol.gsea(path.num, enh_gp = enh_gp)
  ratio <- nrow(pathway.res$enh) / nrow(pathway.res$rna)
  return(ratio)
}

chrom.mod.ratio <- get.ratio(5, df8.all)
 
gsea.ratio <- data.frame(gsea.res.select[, 1])
colnames(gsea.ratio) <- "pathway"

for(i in 1:9){
  ratio <- get.ratio(i, INS.PNET.specSE.enh)
  gsea.ratio[i,9] <- ratio
}

colMeans(gsea.ratio[sapply(gsea.ratio, is.numeric)])

colnames(gsea.ratio)[8] <- "gp2b.ratio"
gsea.ratio[, 2:8] <- round(gsea.ratio[, 2:8], 2) 


pathway.vs.SE <- function(path){
      pathway <- enh.ol.gsea(path, enh_gp = INS.specSE.enh)  
      path.vs.SE <- findOverlaps(toGRanges(pathway$rna), toGRanges(INS.specSE.enh))

      pathway.SE <- cbind(data.frame(chrom.mod.INSspecSE$rna[queryHits(chrom.mod.vs.SE), ]),
                      data.frame(INS.specSE.enh[subjectHits(chrom.mod.vs.SE), ]))

      pathway.SE.sig <- subset(pathway.SE, pathway.SE$FDR < 0.05)
      pathway.sig.cut <- pathway.SE.sig[, c(4,8,10,12,19)]
}

HDM5.vs.SE <- pathway.vs.SE(3)



### get genes from all pathways with high SI 

grn <- rbind(RHO_GTP.all2.pc, RTK.all2.pc, NS.dev.all2.pc, VEGF.all2.pc)

grn.hiSI <- subset(grn, grn$SI >= 4)

write.csv(grn.hiSI[!duplicated(grn.hiSI$gene_name), ]$gene_name, file = "grn_genes.csv")

#### GREAT analysis
library(rGREAT)

chrom.mod.df8.great <- submitGreatJob(chrom.mod.all.enh[, 1:3], species = 'hg38')
avail.ont <- availableOntologies(chrom.mod.df8.great)

df1.great <- submitGreatJob(df1.all[, 1:3], species = 'hg38')
df1.great.enrich <- getEnrichmentTables(df1.great)
df1.great.mf <- df1.great.enrich$`GO Molecular Function`

chrom.mod.great.tb <- getEnrichmentTables(chrom.mod.df8.great)
chrom.mod.mol.func <- chrom.mod.great.tb$`GO Molecular Function`


## SE sig enh
chrom.mod.SE.sig.cut.great <- submitGreatJob(chrom.mod.SE.sig.cut.enh[, 1:3], species = 'hg38')
chrom.mod.SE.tb <- getEnrichmentTables(chrom.mod.SE.sig.cut.great)
