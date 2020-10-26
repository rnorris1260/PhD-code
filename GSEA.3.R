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
TX_biotype <- data.frame(table(rna.res$TX_biotype)[rna.res$TX_biotype])
TX_biotype <- TX_biotype[!duplicated(TX_biotype$Var1), ]

blacklist <- c("nonsense_mediated_decay", "non_stop_decay", "TEC")

rna.res <- na.omit(rna.res)
rna.res.INS_HI.exact <- na.omit(rna.res.INS_HI.exact)

rna.res <- rna.res[rna.res$TX_biotype %ni% blacklist, ]

## select genes and log2FC
all.genes <- rna.res[, c(8,1)]
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
gsea.res.exact <- GSEA(gene_list = ranks, GO_file = "~/pathways/gmt_files/c2.cp.reactome.v7.1.symbols.gmt", pval = 0.05)

## plot
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


####

## select genes from gsea pathways
function(gsea, r){
  genes <- data.frame(gsea[r,8]) 
  colnames(genes) <- "genes"
  rna.res <- tss.ext.qlf.cut[tss.ext.qlf.cut$gene_name %in% genes$genes, ]
  rna.res <- rna.res[order(rna.res$logFC, decreasing = T), ]
  rna.res.genes <- rna.res[!duplicated(rna.res$gene_name), ]
}




## get gene-enhancer associaions

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
  genes.pathway <- rna.res.INS_HI.qlf.canon[rna.res.INS_HI.qlf.canon$gene_name %in% genes$gene_name, 
  return(genes.pathway)
}



#### CREs with mutations
mutations <- read.csv(file = "mutations.csv", header = T)
mutations$chr <- paste0("chr", mutations$chr)
rownames(chrom.mod.LE.mut) <- 1:nrow(chrom.mod.LE.mut)
mut.dist <- distanceToNearest(toGRanges(mutations), toGRanges(rna.res.INS_HI.qlf.canon[, c(9,13,14)]))
mut.dist.df <- data.frame(mut.dist)

mut.dist.df.TSS <- subset(mut.dist.df, mut.dist.df$distance == 0)

mut.TSS <- rna.res[rownames(rna.res) %in% mut.dist.df.TSS$subjectHits, ]
mut.all <- rna.res[rownames(rna.res) %in% mut.dist.df$subjectHits, ]
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
