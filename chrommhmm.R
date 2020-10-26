######## HI H3K27me3 / CHROMHMM

## Liftover to hg38
library(rtracklayer)
#http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
hg19hg38.chain = import.chain("~/tools/hg19ToHg38.over.chain")

hi_hmm <- unlist(liftOver(map, chain))
hi_hmm.df <- data.frame(hi_hmm)

hg18.chain <- import.chain("~/tools/hg18ToHg38.over.chain")

setwd("~/data/H3K27me3/")
HI.H3K27me3 <- read.bed("HI_H3K27me3.cps.bed")

HI.H3K27me3 <- regioneR::toGRanges(HI.H3K27me3) 
HI.H3K27me3 <- unlist(liftOver(HI.H3K27me3, hg19hg38.chain))



####### prepare enh data
#ins <- regioneR::toGRanges(enh)

#chips.res.dv <- chips.res.merge[, c(1:4,6,7,9,10,11)]


df1.all <- subset(chip.fc.fdr.all, chip.fc.fdr.all$type == "DI-2")
df1b[, 8] <- "1"
df2.all <- subset(chip.fc.fdr.all, chip.fc.fdr.all$type == "DI-HI")
df2b[, 8] <- "2"
df3.all <- subset(chip.fc.fdr.all, chip.fc.fdr.all$type == "DPNET")
df3b[, 8] <- "3"
df4.all <- subset(chip.fc.fdr.all, chip.fc.fdr.all$type == "D2")
df4b[, 8] <- "4"
df5.all <- subset(chip.fc.fdr.all, chip.fc.fdr.all$type == "DHI")
df5b[, 8] <- "5"
df6.all <- subset(chip.fc.fdr.all, chip.fc.fdr.all$type == "DI-PNET")
df6b[, 8] <- "6"

df7.all <- rbind(df4.all, df5.all)
df7b[, 8] <- "7"
df8.all <- rbind(df1.all, df2.all)
df8b[, 8] <- "8"
df9 <- subset(chip.fc.fdr, chip.fc.fdr$type == "non.sig")
df9b[, 8] <- "9"


# Find overlaps between hi_hmm and enhancers
ac.vs.chromhmm <- function(enh){
      enh.df <- data.frame(enh$chr, enh$start, enh$end, enh$volcano.gp)
      colnames(enh.df) <- c("chr", "start", "end", "volcano.gp")
      rownames(enh.df) <- c(1:nrow(enh.df))
      enh.gr <- regioneR::toGRanges(enh.df)
      ols <- findOverlaps(enh.gr, hi_hmm)

      df <- cbind(data.frame(mcols(enh.gr))[queryHits(ols),],
            data.frame(mcols(hi_hmm))[subjectHits(ols),])
      df <- df[, c(1,4)]
      colnames(df) <- c("volcano.gp", "ChromHMM_annotation")
      df <- df

      df.nohits <- enh[rownames(enh.df) %ni% queryHits(ols), ]
      df.nohits <- data.frame(df.nohits$volcano.gp)
      colnames(df.nohits) <- "volcano.gp"
      df.nohits$ChromHMM_annotation <- "No hit"

      df.all <- rbind(df, df.nohits)

}




##### loop

dfs <- list(df1.all, df2.all, df3.all, df4.all, df5.all, df6.all, 
            df7.all, df8.all, PNET_HI_DA, df9)


chromm.hmm.dv <- data.frame()

for(i in dfs){
    df <- ac.vs.chromhmm(i)
    chromm.hmm.dv <- rbind(chromm.hmm.dv, df)
}


chromm.hmm.dv$ChromHMM_annotation <- factor(chromm.hmm.dv$ChromHMM_annotation, 
                                            levels=c("Active promoter",
                                                     "Weak promoter",
                                                     "Poised promoter",
                                                     "Transcribed",
                                                     "Strong enhancer",
                                                     "Weak enhancer",
                                                     "Insulator",
                                                     "Repressed",
                                                     "No hit"))

#chip.res.me3 <- ac.vs.me3(chips.res.merge)
#colnames(chromm.hmm.res9)[7] <- "type"
#cols <- map_col


cols <- c("light blue", "pink", "chartreuse2",
                 "yellow", "dark orange", "purple",
                  "red", "grey", "gray20")
names(cols) <- c("Insulator", "Weak promoter", "Transcribed", "Weak enhancer", "Strong enhancer", "Poised promoter",
                     "Active promoter", "Repressed", "No hit")
#colScale <- scale_colour_manual(name = "ChromHMM_Annotation", values = myColors)

#col.list <- split(col.df, seq(nrow(col.df)))

library(ggplot2)

setwd("~/plots/thesis/")

png(filename='chromHMM.png', height = 700, width=700)

ggplot(chromm.hmm.dv,
       aes(volcano.gp, ..count..)) +
  geom_bar(aes(fill=factor(ChromHMM_annotation)), position="fill", color="black") +
  scale_fill_manual(values=cols,
                    name=expression("HI ChromHMM")) +
  scale_y_continuous(name="Percentage (%)",
                     labels=function(x) x*100) +
  xlab("H3K27ac type") +
  cowplot::theme_cowplot() +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()


##### get regions in section 8 of double volcano that are repressed in human islets
df8.rep <- subset(chromm.hmm.res, chromm.hmm.res$type == 8 & chromm.hmm.res$ChromHMM_Annotation == "Repressed")
df8.rep.bed <- chips.res.merge[chips.res.merge$id %in% df8.rep$id, ]
df8.rep.gr <- regioneR::toGRanges(df8.rep.bed[, 1:3])



    

      
