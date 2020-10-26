
library(GenomicRanges)
library(regioneR)


#### investigate individual hits from dv group
setwd("~/data/H3K27me3/original_procesed/")
HI.H3K27me3.1 <- read.bed("BI.Pancreatic_Islets.H3K27me3.pancreatic_islets_normal_3_27_09.hg38_peaks.broadPeak")
HI.H3K27me3.1 <- HI.H3K27me3.1[HI.H3K27me3.1$chr %in% chrs$chr, ]
HI.H3K27me3.1 <- toGRanges(HI.H3K27me3.1[, 1:3])

HI.H3K27me3.2 <- read.bed("UCSF-UBC.Pancreatic_Islets.H3K27me3.ZGI_213.hg38_peaks.broadPeak")
HI.H3K27me3.2 <- HI.H3K27me3.2[HI.H3K27me3.2$chr %in% chrs$chr, ]
HI.H3K27me3.2 <- toGRanges(HI.H3K27me3.2[, 1:3])

HI.H3K27me3 <- regioneR::joinRegions(c(HI.H3K27me3.1, HI.H3K27me3.2))


get.me3 <- function(enh){
        enh.gr <- toGRanges(enh)
        ols.me3 <- findOverlaps(enh.gr, HI.H3K27me3)

        enh$H3K27me3 <- "No hit"
        me3.hits <- unique(queryHits(ols.me3))
        enh$H3K27me3[me3.hits] <- "HI-H3K27me3"

return(enh)
}
  
#df1.all, df2.all
dv.groups <- rbind(df1.all, df2.all, df3.all, df4.all, df5.all, df6.all, df7.all, df8.all, df9)
dfs <- list(df1.all, df2.all, df3.all, df4.all, df5.all)


dv.groups.me3 <- get.me3(dv.groups)
colnames(dv.groups.me3)[8] <- "volcano_group"

PNET_HI_DA.me3 <- get.me3(PNET_HI_DA)
PNETnb.DA.me3 <- get.me3(PNETnb.DA)

### get % overlap for each dv group

#cps.enhtype.me3 <- get.me3(cps.enhtype)
#cps.enhtype.beta.me3 <- get.me3(cps.enhtype.beta)


#dv.gps.me3.hits <- subset(dv.groups.me3, dv.groups.me3$H3K27me3 == "HI-H3K27me3")


me3.df.new <- data.frame()

for(i in unique(dv.groups.me3$volcano.gp)){
  data <- subset(dv.groups.me3, dv.groups.me3$volcano.gp == i)
  total <- nrow(data)
  me3 <- nrow(subset(data, data$H3K27me3 == "HI-H3K27me3"))
  me3 <- (100/total) * me3
  nohit <- nrow(subset(data, data$H3K27me3 == "No hit"))
  nohit <- (100/total * nohit)
  df <- data.frame(i, me3, nohit)
  me3.df.new <- rbind(me3.df.new, df)
}

me3.indiv <- function(data){
    total <- nrow(data)
    me3 <- nrow(subset(data, data$H3K27me3 == "HI-H3K27me3"))
    me3 <- (100/total) * me3
    nohit <- nrow(subset(data, data$H3K27me3 == "No hit"))
    nohit <- (100/total * nohit)
    df <- data.frame(me3, nohit)
}

PNET_HI_DA.me3.res <- me3.indiv(PNET_HI_DA.me3)
PNETnb.DA.me3.res <- me3.indiv(PNETnb.DA.me3)


library(reshape)
colnames(me3.df)[2] <- "H3K27me3"

me3.plot <- melt(me3.df, id.vars = "volcano_group") 

colnames(me3.plot) <- c("volcano_group", "group", "value")

setwd("~/plots/thesis/")

png(filename='ac_vs_me3.png', height = 700, width=700)
ggplot(me3.plot[1:7, ], aes(x=factor(volcano_group), y=value)) +
     geom_bar(stat = "identity", fill = "lightcoral") +
      facet_wrap(~group, ncol = 1, scales = "free") +
        xlab("Volcano group") +
        ylab("Percentage overlap") +
     theme(legend.position="top", legend.text = element_text(size = 20), 
           axis.text = element_text(size = 20), axis.title = element_text(size = 20))
dev.off()



### permutation test
genome <- getGenome("hg38")

perm.res2 <- data.frame()

for(i in dfs){
  df.perm <- toGRanges(i)
  pt <- overlapPermTest(A=df.perm, B=HI.H3K27me3, ntimes=1000, genome = genome,
                        alternative = "greater", verbose = TRUE)
  volc.gp <- unique(i$type)
  pvalue <- pt$numOverlaps$pval
  zscore <- pt$numOverlaps$zscore
  res <- data.frame(volc.gp, pvalue, zscore)
  perm.res2 <- rbind(perm.res2, res)
}


PNETnb.DA.pt <- overlapPermTest(A=toGRanges(PNETnb.DA), B=HI.H3K27me3, ntimes=1000, genome = genome,
                          alternative = "greater", verbose = TRUE)

perm.res2$FDR <- p.adjust(perm.res2$pvalue, method = "bonferroni", n = 5)



ggplot(data.frame(dv.groups.me3),
       aes(volcano_group)) +
  geom_bar(aes(H3K27me3), position="fill", color="black") +
  xlab("H3K27ac type") +
  scale_y_continuous(name="Percentage (%)",
                     labels=function(x) x*100) +
  cowplot::theme_cowplot() +
  theme(legend.position="top")



dv.group.1 <- subset(dv.groups.me3, dv.groups.me3$volcano_group == "1" &
                       dv.groups.me3$H3K27me3 == "HI-H3K27me3")


ggplot(data.frame(cps.enhtype.beta.me3),
       aes(type)) +
  geom_bar(aes(fill=H3K27me3), position="fill", color="black") +
  xlab("H3K27ac type") +
  scale_y_continuous(name="Percentage (%)",
                     labels=function(x) x*100) +
  coord_flip() + 
  cowplot::theme_cowplot() +
  theme(legend.position="top")




ggplot(cps.enhtype.me3,
       aes(type, ..count..)) +
  geom_bar(aes(fill=H3K27me3), position="fill", color="black") +
  scale_fill_manual(values=c("red", "blue"),
                    name=expression("H3K27me3")) +
  scale_y_continuous(name="Percentage (%)",
                     labels=function(x) x*100) +
  xlab("H3K27ac type") +
  cowplot::theme_cowplot()





ac.vs.me3 <- function(enh){
  enh.gr <- regioneR::toGRanges(enh)
  ols.me3 <- findOverlaps(enh.gr, HI.H3K27me3)
  
  enh$H3K27me3 <- "No hit"
  enh$H3K27me3[queryHits(ols.me3)] <- "H3K27me3"
  
  me3.res <- enh
}






k27me3.res <- data.frame()

for(i in dfs){
  df.me3 <- ac.vs.me3(i)
  k27me3.res <- rbind(k27me3.res, df.me3)
}


k27me3.res$volcano_group <- as.character(k27me3.res$volcano_group)
k27me3.res$volcano_group <- factor(k27me3.res$volcano_group, level=c("1", "4", "2", "5", "6", "3", "8", "7", "9-.01", "9-.05"))


k27me3.res9.v2 %>% 
  dplyr::mutate(type = factor(type, levels = c("1", "4", "2", "5", "6", "3", "8", "7", "9-.01", "9-.05"))) %>% 
  k27me3.res9.v2







###### convert section dfs to 3 column bed files
df1.bed <- df1[, c(1:3,10)]
df2.bed <- df2[, c(1:3,10)]
df8.bed <- df8[, c(1:3,10)]
df1.bed[, 4] <- "group-1"
df2.bed[, 4] <- "group-2"
df8.bed[, 4] <- "group-8"

setwd("~/data/ChIP/bed_files/")

INS.cps <- read.bed("INS12.cps.bed")
HI6.cps <- read.bed("HI6.cps.bed")
PNET.cps <- read.bed("PNET.cps.bed")

INS.cps[, 4] <- "INS.cps"
HI6.cps[, 4] <- "HI.cps"
PNET.cps[, 4] <- "PNET.cps"







dfs.cps <- list(df1.bed, df2.bed, df8.bed, INS.cps, HI6.cps, PNET.cps)
dfs.cps <- list(INS.cps, HI6.cps, PNET.cps)

k27me3.res.cps <- data.frame()
k27me3.res.cps2 <- data.frame()

for(i in dfs.cps){
  df.me3 <- ac.vs.me3.bed(i)
  k27me3.res.cps <- rbind(k27me3.res.cps, df.me3)
}

colnames(k27me3.res.cps)[c(5,6)] <- c("H3K27me3", "enh_type")


ggplot(k27me3.res.cps[, c(5,6)],
       aes(enh_type, ..count..)) +
  geom_bar(aes(fill=H3K27me3), position="fill", color="black") +
  scale_fill_manual(values=c("red", "blue"),
                    name=expression("H3K27me3")) +
  scale_y_continuous(name="Percentage (%)",
                     labels=function(x) x*100) +
  xlab("H3K27ac type") +
  cowplot::theme_cowplot()





Ttags.intersect <- read.bed("Ttags.INS_HI_PNET.intersect.bed")
Ttags.intersect <- Ttags.intersect[, c(1:3,8)]
Ttags.intersect$region <- paste(Ttags.intersect$chr, Ttags.intersect$start, Ttags.intersect$end)

colnames(Ttags.intersect)[4] <- "cps"

Ttags.intersect <-Ttags.intersect[!duplicated(Ttags.intersect[, 1:5]), ]

hits2 <- data.frame()

for(i in unique(Ttags.intersect2$region)){
  
  hit <- subset(Ttags.intersect2, Ttags.intersect2$region == i)
  
  if (nrow(hit) == 3) {hit$type <- "INS/HI/PNET"}
  
  if (nrow(hit) == 1 && hit$cps == "1") {hit$type <- "INS"}
  if (nrow(hit) == 1 && hit$cps == "2") {hit$type <- "HI"}
  if (nrow(hit) == 1 && hit$cps == "3") {hit$type <- "PNET"}
  
  if (nrow(hit) == 2 && hit[1, 4] == "1" && hit[2, 4] == "2" | 
      nrow(hit) == 2 && hit[1, 4] == "2" && hit[2, 4] == "1") {hit$type <- "INS/HI"}
  
  if (nrow(hit) == 2 && hit[1, 4] == "1" && hit[2, 4] == "3" | 
      nrow(hit) == 2 && hit[1, 4] == "3" && hit[2, 4] == "1") {hit$type <- "INS/PNET"}
  
  if (nrow(hit) == 2 && hit[1, 4] == "2" && hit[2, 4] == "3" | 
      nrow(hit) == 2 && hit[1, 4] == "3" && hit[2, 4] == "2") {hit$type <- "HI/PNET"}
  
  hit <- hit[!duplicated(hit$region), ]
  hits2 <- rbind(hits2, hit)  
}


hit <- Ttags.intersect[1:2, ]
cps <- list(hit$cps)

cps.enhtype <- rbind(hits, hits2)
cps.enhtype <- cps.enhtype[!duplicated(cps.enhtype$region), ]
cps.enhtype <- cps.enhtype[, -4]

rownames(cps.enhtype) <- c(1:nrow(cps.enhtype))

cps.ols <- findOverlaps()

ins$H3K27me3 <- "No"
ins$H3K27me3[queryHits(ols)] <- "H3K27me3"



k27me3.cps.combos <- data.frame()

for(i in unique(cps.enhtype$type)){
  df.cps <- subset(cps.enhtype, cps.enhtype$type == i)
  df.me3 <- ac.vs.me3(df.cps)
  k27me3.cps.combos <- rbind(k27me3.cps.combos, df.me3)
}


ggplot(k27me3.cps.combos,
       aes(type, ..count..)) +
  geom_bar(aes(fill=H3K27me3), position="fill", color="black") +
  scale_fill_manual(values=c("red", "blue"),
                    name=expression("H3K27me3")) +
  scale_y_continuous(name="Percentage (%)",
                     labels=function(x) x*100) +
  xlab("H3K27ac cps type") +
  theme(axis.text.x=element_text(angle=90, hjust=1))
cowplot::theme_cowplot()



##### PNETs...

cps.enhtype.PNET <- subset(cps.enhtype, cps.enhtype$type == "INS/PNET" | cps.enhtype$type == "HI/PNET" |
                             cps.enhtype$type == "INS/HI/PNET" | cps.enhtype$type == "PNET")

write.bed(cps.enhtype.PNET[, 1:3], file = "cps_enhtype_PNET.bed")


PNETs_type.intersect <- read.bed("PNET_types.int.bed")  

### 1 = alpha , 2 = beta, 3 = DP

PNETs_type.intersect <- PNETs_type.intersect[!duplicated(PNETs_type.intersect[, 1:4]), ]

PNET_alpha.cps <- read.bed("PNET_alpha.cps.bed")
PNET_beta.cps <- read.bed("PNET_beta.cps.bed")
PNET_DP.cps <- read.bed("PNET_DP.cps.bed")
INS.cps <- read.bed("INS12.cps.bed")
HI.cps <- read.bed("HI6.cps.bed")

HI.cps$type <- "HI"



cps.list <- list(PNET_alpha.cps, PNET_beta.cps, PNET_DP.cps, INS.cps, HI.cps)

##### add in DARE-INS and all chip results without DARE-INS

DARE_INS.bed$type <- "DARE-INS"
topTags.chip.minDI$type <- "NOT-DI"
cps.list <- list(PNET_alpha.cps, PNET_beta.cps, PNET_DP.cps, INS.cps, HI.cps, DARE_INS.bed, topTags.chip.minDI)


k27me3.INS.HI.PNETbeta.DI <- data.frame()

for(i in cps.list){
  df.me3 <- ac.vs.me3.bed(i)
  k27me3.INS.HI.PNETbeta.DI <- rbind(k27me3.INS.HI.PNETbeta.DI, df.me3)
}

colnames(k27me3.INS.HI.PNETbeta)[c(5,6)] <- c("H3K27me3", "enh_type")







for(i in unique(cps.enhtype$type)){
  df.cps <- subset(cps.enhtype, cps.enhtype$type == i)
  df.me3 <- ac.vs.me3.bed(df.cps)
  k27me3.INS.HI.PNETbeta <- rbind(k27me3.INS.HI.PNETbeta, df.me3)
}

colnames(k27me3.INS.HI.PNETbeta.DI)[c(5,6)] <- c("H3K27me3", "enh_type")

ggplot(k27me3.INS.HI.PNETbeta.DI[, c(5,6)],
       aes(volcano_group, ..count..)) +
  geom_bar(aes(fill=H3K27me3), position="fill", color="black") +
  scale_fill_manual(values=c("red", "blue"),
                    name=expression("H3K27me3")) +
  scale_y_continuous(name="Percentage (%)",
                     labels=function(x) x*100) +
  xlab("H3K27ac type") +
  theme(axis.text.x=element_text(angle=90, hjust=1))
cowplot::theme_cowplot()


PNETs_type.intersect$region <- paste(PNETs_type.intersect$chr, PNETs_type.intersect$start, 
                                     PNETs_type.intersect$end)
PNETs_type.intersect[, c(4,5)] <- table(PNETs_type.intersect$region)[PNETs_type.intersect$region]

colnames(PNETs_type.intersect)[5] <- "freq"

PNETs_type.beta.spec <- subset(PNETs_type.intersect, PNETs_type.intersect$freq == 1 &
                                 PNETs_type.intersect$PNETs_type == "2")


cps.enhtype.PNET.gr <- toGRanges(cps.enhtype.PNET)
PNET_beta.cps.gr <- toGRanges(PNET_beta.cps)

ols <- findOverlaps(cps.enhtype.PNET.gr, PNET_beta.cps.gr)
cps.enhtype.beta <- cps.enhtype.PNET[queryHits(ols), ] 

cps.enhtype.INS.HI <- subset(cps.enhtype, cps.enhtype$type == "INS" | cps.enhtype$type == "HI" |
                               cps.enhtype$type == "INS/HI")


cps.enhtype <- rbind(cps.enhtype.INS.HI, cps.enhtype.beta)




####### so..... Insulinoma (but not PNET) development derives from the activation poised HI enhancers

#### so let's select all regions from groups 8 (1 + 2) that overlap with H3K27me3, find corresponding nfr regions and then get motifs

DI.g8.me3 <- subset(k27me3.res9.v2, k27me3.res9.v2$type == "8" & k27me3.res9.v2$H3K27me3 == "HI-H3K27me3")







######### super-enhancers

## all chip res regions
SE.1 <- chip.res.todo[, 1:3]
SE.1$type <- "SE.all"


## regions that overlap an HI SE
SE.2 <- subset(chip.res.todo, chip.res.todo$HI_SE != "none")
SE.2 <- SE.2[, 1:3]
SE.2$type <- "all.HI"


## regions that do not overlap an HI SE
SE.3 <- subset(chip.res.todo, chip.res.todo$HI_SE == "none")
SE.3 <- SE.3[, 1:3]
SE.3$type <- "no.HI"


## regions that overlap an INS SE
SE.4 <- subset(chip.res.todo, chip.res.todo$INS_SE != "none")
SE.4 <- SE.4[, 1:3]
SE.4$type <- "all.INS"

## regions that do not overlap an INS SE
SE.5 <- subset(chip.res.todo, chip.res.todo$INS_SE == "none")
SE.5 <- SE.5[, 1:3]
SE.5$type <- "no.INS"

## INS-specific
SE.6 <- subset(chip.res.todo, chip.res.todo$INS_SE != "none" &
                 chip.res.todo$HI_SE == "none")
SE.6 <- SE.6[, 1:3]
SE.6$type <- "INS-specific"

## HI-specific
SE.7 <- subset(chip.res.todo, chip.res.todo$HI_SE != "none" &
                 chip.res.todo$INS_SE == "none")
SE.7 <- SE.7[, 1:3]
SE.7$type <- "HI-specific"


SE.dfs <- list(SE.1, SE.2, SE.3, SE.4, SE.5, SE.6, SE.7)


k27me3.SEs <- data.frame()

for(i in SE.dfs){
  df.me3 <- ac.vs.me3.bed(i)
  k27me3.SEs <- rbind(k27me3.SEs, df.me3)
}

colnames(k27me3.SEs)[5] <- c("H3K27me3")


ggplot(k27me3.SEs[, c(5,6)],
       aes(type, ..count..)) +
  geom_bar(aes(fill=H3K27me3), position="fill", color="black") +
  scale_fill_manual(values=c("red", "blue"),
                    name=expression("H3K27me3")) +
  scale_y_continuous(name="Percentage (%)",
                     labels=function(x) x*100) +
  xlab("H3K27ac type") +
  theme(axis.text.x=element_text(angle=90, hjust=1))
cowplot::theme_cowplot()


