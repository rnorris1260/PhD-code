#########

######## select enhancer regions
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
DI.nfr.ATAC <- ATACvsDI(ATAC.bed = "INS_nfr_ATAC.cps.bed", DI.gr = DARE_INS.gr)










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

  








  
