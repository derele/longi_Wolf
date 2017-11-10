library(MultiAmplicon)
library(dada2)
library(parallel)

FILTER <- FALSE
newMA <- FALSE


a.files <- list.files(path="/SAN/Ines_longitudinal/NGS3/results/Unaligned_Mismatch0/fastq files/",
                      pattern=".fastq.gz", full.names=TRUE)

Ffq.file <- a.files[grepl("R1", a.files)]
Rfq.file <- a.files[grepl("R2", a.files)]

## asses the Quality
## plotQualityProfile(Rfq.file[[20]])

samples <- gsub("/SAN/Ines_longitudinal/NGS3/results/Unaligned_Mismatch0/fastq files//(.*?S\\d+).*?\\.fastq\\.gz", "\\1", Ffq.file)

filt_path <- "/SAN/Ines_longitudinal/NGS3/results/filt_fastq/"
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(samples, "_F_filt.fastq.gz"))
names(filtFs) <- samples

filtRs <- file.path(filt_path, paste0(samples, "_R_filt.fastq.gz"))
names(filtRs) <- samples

if (FILTER){
    ## Filter # only run when new filtered data is needed
    mclapply(seq_along(Ffq.file), function(i) {
        fastqPairedFilter(c(Ffq.file[i], Rfq.file[i]), c(filtFs[i], filtRs[i]),
                          truncLen=c(170,170), 
                          maxN=0, maxEE=2, truncQ=2, 
                          compress=TRUE, verbose=TRUE)
    }, mc.cores=20)
}



file.ptable <- "/SAN/Ines_longitudinal/NGS3/results/Unaligned_Mismatch0/primers_comb.csv"
ptable <- read.csv(file.ptable, sep=",", header=FALSE)

primerF <- as.character(ptable[,2])
primerR <- as.character(ptable[,3])

names(primerF) <- paste0(ptable[,1], "_F")
names(primerR) <- paste0(ptable[,1], "_R")

rownames(ptable) <- ptable[,1]

files <- PairedReadFileSet(filtFs, filtRs)
primers <- PrimerPairsSet(primerF, primerR)

if (newMA){
    MA <- MultiAmplicon(primers, files)
    MA1 <- sortAmplicons(MA)    
    pdf("figures/primers_MA_sorted.pdf", 
        width=25, height=15, onefile=FALSE)
    cluster <- plot_Amplicon_numbers(rawCounts(MA1))
    dev.off()
    MA2 <- derepMulti(MA1)
    MA3 <- dadaMulti(MA2, err=NULL, selfConsist=TRUE,
                     multithread=TRUE)
    MA4 <- MultiAmplicon:::mergeMulti(MA3, justConcatenate=TRUE)
    MA5 <- MultiAmplicon:::sequenceTableMulti(MA4)
    MA6 <- MultiAmplicon:::noChimeMulti(MA5, mc.cores=20)
    names(MA6@sequenceTableNoChime) <- rownames(MA6)
    STnoC <- MA6@sequenceTableNoChime
    save(STnoC, file="/SAN/Ines_longitudinal/Seqtable.Rdata")
} else {load(file="/SAN/Ines_longitudinal/Seqtable.Rdata")} ## -> STnoC

