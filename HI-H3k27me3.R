
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









