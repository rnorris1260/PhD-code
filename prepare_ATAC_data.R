#########

######## select enhancer regions
DARE_INS <- DARE_INS[order(DARE_INS$FDR), ]
load("DI.HI.rda")


ATACvsDI <- function(ATAC.bed, DI.gr){
## read bed
ATAC.cps <- read.bed(ATAC.bed)
## remove xs chrs
ATAC.cps <- ATAC.cps[ATAC.cps$chr %in% chrs$chr, ]
## convert to gr
ATAC.gr <- regioneR::toGRanges(ATAC.cps)
## find overlaps
DIvsATAC <- findOverlaps(ATAC.gr, DI.gr)
## create df
DI.open <- cbind(data.frame(ATAC.gr[queryHits(DIvsATAC),]),
                    data.frame(DI.gr[subjectHits(DIvsATAC),]))
## return df with uniwue NFR hits
DI.open <- DI.open[!duplicated(DI.open[, 1:3]), ]
}


setwd("~/data/ATAC/")
DI.nfr <- ATACvsDI(ATAC.bed = "INS_nfr_cps.bed", DI.gr = DARE_INS.gr)
DI.ATAC <- ATACvsDI(ATAC.bed = "INS_beta_HI_PC_ATAC.cps.bed", DI.gr = DARE_INS.gr)
DI.HI.opan.bed <- DI.HI.opan[!duplicated(DI.HI.opan[, 1:3]), ]
DI.nfr.ATAC <- ATACvsDI(ATAC.bed = "INS_nfr_ATAC.cps.bed", DI.gr = DARE_INS.gr)

gp1_all.nfr <- ATACvsDI(ATAC.bed = "INS_nfr_ATAC.cps.bed", DI.gr = toGRanges(df1.all))
gp1_all.nfr <- gp1_all.nfr[!duplicated(gp1_all.nfr[, 1:3]), ]

gp8_all.nfr <- ATACvsDI(ATAC.bed = "INS_nfr_ATAC.cps.bed", DI.gr = toGRanges(df8.all))
gp8_all.nfr <- gp8_all.nfr[!duplicated(gp8_all.nfr[, 1:3]), ]


DI.g8.me3.gr <- toGRanges(DI.g8.me3[, 1:3])
DI.g8me3.nfr <- ATACvsDI(ATAC.bed = "INS_nfr_ATAC.cps.bed", DI.gr = DI.g8.me3.gr)
gp8.nfr <- ATACvsDI(ATAC.bed = "INS_nfr_ATAC.cps.bed", DI.gr = toGRanges(df8b))

dv2.in.gsea.nfr <- ATACvsDI(ATAC.bed = "INS_nfr_ATAC.cps.bed", DI.gr = toGRanges(dv2.in.gsea.enh))

dv2.nfr <- ATACvsDI(ATAC.bed = "INS_nfr_ATAC.cps.bed", DI.gr = toGRanges(df2b))




#### df 1 with HI logFC > 2
df1.fc2.1.nfr <- ATACvsDI(ATAC.bed = "INS_nfr_ATAC.cps.bed", DI.gr = dv.gp1.fc2)

## DI.PNET
DI.PNET <- subset(chips.res.merge, chips.res.merge$enh_type == "DI.PNET")
DI.PNET.gr <- regioneR::toGRanges(DI.PNET[, 1:3])    
DI.PNET.nfr.ATAC <- ATACvsDI(ATAC.bed = "INS_nfr_ATAC.cps.bed", DI.gr = DI.PNET.gr)



setwd("~/motifs/")
write.bed(gp8_all.nfr[, 1:3], file = "gp8_all_nfr.bed")



INSbetaHI.ATAC <- read.bed("INS_beta_HI_ATAC.cps.bed")
INSbetaHI.ATAC <- INSbetaHI.ATAC[INSbetaHI.ATAC$chr %in% chrs$chr, ]
INSbetaHI.ATAC.gr <- regioneR::toGRanges(INSbetaHI.ATAC)

INS_nfr.ATAC.gr <- regioneR::toGRanges(INS_nfr.cps.bed) 

nfr.vs.ATAC <- findOverlaps(INSbetaHI.ATAC.gr, INS_nfr.ATAC.gr)
nfr.vs.ATAC.df <- data.frame(nfr.vs.ATAC)


nfr.vs.ATAC.dist <- data.frame(distanceToNearest(INSbetaHI.ATAC.gr, INS_nfr.ATAC.gr))



##### get enhancer res with sharing index
INS_HI.chip.res <- INS.HI.chip.res[, c(1:3,22:26)]
INS_REs <- read.bed("INS_REs.bed")
colnames(INS_REs)[c(4,5)] <- c("enh_id", "INS.sharing.index")
INS.chip.res <- dplyr::left_join(INS_HI.chip.res, INS_REs)


#### select DARE_INS with higher sharing index - but not 12 as that will select for contaminants
DI.sigSI <- subset(INS.chip.res, 
                   INS.chip.res$logFC > 1 &
                    INS.chip.res$FDR < 0.05 &
                   INS.chip.res$INS.sharing.index >2 &
                     INS.chip.res$INS.sharing.index < 12)
## get nfr
DI.sigSI.gr <- regioneR::toGRanges(DI.sigSI[, 1:3])
DIsigSI.nfr.ATAC <- ATACvsDI(ATAC.bed = "INS_nfr_ATAC.cps.bed", DI.gr = DI.sigSI.gr)

write.bed(DIsigSI.nfr.ATAC[, 1:3], file = "DIsigSI.nfr.ATAC.bed")


#### select super-enhancers

INS_SE.enh <- subset(chip.res.todo, chip.res.todo$INS_SE != "none" &
                       chip.res.todo$HI_SE == "none")

INS_SE.enh.gr <- regioneR::toGRanges(INS_SE.enh[, 1:3])
INS_SE.nfr.ATAC <- ATACvsDI(ATAC.bed = "INS_nfr_ATAC.cps.bed", DI.gr = INS_SE.enh.gr)
write.bed(INS_SE.nfr.ATAC, file = "INS_SE_nfr_ATAC.bed")

  

###### are there more overlaps (H3K27ac vs ATAC) for super-enhancers vs regular enhancers???

### select enhancers present in at least one insulinoma
INS_REs <- subset(chip.res.todo, chip.res.todo$INS_SI <=1)

### separate df into SE and non-SE  
INS_RE_noSE <- subset(INS_REs, INS_REs$INS_SE == "none")
INS_RE_SE <- subset(INS_REs, INS_REs$INS_SE != "none")

ATACvsenh <- function(enh_list){
  enh.gr <- regioneR::toGRanges(enh_list)
  enhvsATAC <- findOverlaps(INS_PC_ATAC.hg38.gr, enh.gr)
  ## create df
  enh.open <- cbind(data.frame(INS_PC_ATAC.hg38.gr[queryHits(enhvsATAC),]),
                      data.frame(enh.gr[subjectHits(enhvsATAC),]))
  ## get number of enhancers with corresponding ATAC peak
  #return(nrow(unique(data.frame(enh.open[, 6:8]))))
}

INS_RE_noSE.open <- ATACvsenh(INS_RE_noSE)



####### HOM|ER nfr
convert.nfr <- function(txt){ 
      nfr <- read.table(txt, header=F, skip = 34)
      nfr <- nfr[, 2:4]
      colnames(nfr) <- c("chr", "start", "end")
      nfr <- nfr[order(nfr$chr, nfr$start), ]
      name <- gsub("homerpeaks.txt", "nfr.bed", txt)
      write.bed(nfr, file = paste0(nfr, "nfr.bed"))
}


files <- list.files(path = "~/data/ChIP/Motifs/homer_nfr_bed", 
                               pattern = "*.homerpeaks.txt$", full.names = T)
  
files <- c("NET10_INS.homerpeaks.txt", "NET11_INS.homerpeaks.txt", "NET14_INS.homerpeaks.txt", "NET16_INS.homerpeaks.txt",
           "NET17_INS.homerpeaks.txt", "NET20_INS.homerpeaks.txt", "NET21_INS.homerpeaks.txt", "NET25_INS.homerpeaks.txt",
           "NET26_INS.homerpeaks.txt", "NET29INS.homerpeaks.txt", "NET30_INS.homerpeaks.txt", "NET38_INS.homerpeaks.txt")

  
for(i in files){
  convert.nfr(i)
}

  
nfr <- read.table("NET10_INS.homerpeaks.txt", header=F, skip = 34)

quote="\""




INS_nfr.sort.bed <- read.bed("INS_nfr_sort.bed")
INS_nfr.cps.bed <- read.bed("INS_nfr_cps.bed")



#####################

# Get open regions for section 8 of double volcano plot

df8.rep.nfr.ATAC <- ATACvsDI(ATAC.bed = "INS_nfr_ATAC.cps.bed", DI.gr = df8.rep.gr)

write.bed(df8.rep.nfr.ATAC[, 1:3], file = "df8_rep_open.bed")
  