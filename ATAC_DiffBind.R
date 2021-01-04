# Load in libraries
library(DESeq2)
library(BiocParallel)
library(DiffBind)
library(GenomicRanges)

# Set number of processors
register(BPPARAM = MulticoreParam(4))

# Read in data and establish contrast
ATAC <- dba(sampleSheet = "sample_sheet.csv", peakFormat="narrow")
ATAC_count <- dba.count(ATAC,summits=250,bParallel = TRUE)
ATAC_count <- dba.contrast(ATAC_count, categories=c(DBA_TREATMENT,DBA_CONDITION))

# Run DiffBind analysis with edgeR
ATAC_count <-dba.analyze(ATAC_count,method=DBA_EDGER,bParallel = TRUE)

# EB vs. Veh comparison 
ATAC_pairwise_treatment <- dba.report(ATAC_count,contrast = 1,fold=1,th=0.05,method = DBA_EDGER)
write.table(ATAC_pairwise_treatment,"ATAC_EDGER_treatment.txt",quote = FALSE,row.names = FALSE,col.names=FALSE,sep = '\t')

# Male vs. Female comparison
ATAC_pairwise_sex <- dba.report(ATAC_count,contrast = 2,fold=1,th=0.05,method = DBA_EDGER)
write.table(ATAC_pairwise_sex,"ATAC_EDGER_sex.txt",quote = FALSE,row.names = FALSE,col.names=FALSE,sep = '\t')


ATAC_pairwise_female <- dba.report(ATAC_count,contrast = 7,th=0.01,method = DBA_EDGER,bFlip = TRUE,bUsePval=TRUE)
write.table(ATAC_pairwise_male[ATAC_pairwise_male$Fold>0,],"/Users/gegenhu/Desktop/ATAC_EDGER_pairwise_male_pval_0.01.bed",quote = FALSE,col.names=FALSE,row.names = FALSE,sep = '\t')
write.table(ATAC_pairwise_female[ATAC_pairwise_female$Fold>0,],"/Users/gegenhu/Desktop/ATAC_EDGER_pairwise_female_pval_0.01.bed",quote = FALSE,col.names=FALSE,row.names = FALSE,sep = '\t')


ATAC_MA_intervals <- dba.report(ATAC_count,contrast = 1,th=1,method = DBA_EDGER,bFlip = TRUE)
write.table(ATAC_MA_intervals,"/Users/gegenhu/Documents/TOLLKUHN_LAB/ATAC_seq/ATAC_EDGER_MA_intervals.bed",quote = FALSE,row.names = FALSE,sep = '\t')



ATAC_mEBvmCO <- dba.report(ATAC_count,contrast = 4,th=0.05,fold=1,method = DBA_EDGER,bFlip = TRUE)
write.table(ATAC_mEBvmCO[ATAC_mEBvmCO$Fold>0,],"/Volumes/tollkuhn_hpc_norepl_home/gegenhu/Adult_ATAC/ATAC_EDGER_mEBvmCO.bed",quote = FALSE,row.names = FALSE,col.names=FALSE,sep = '\t')
