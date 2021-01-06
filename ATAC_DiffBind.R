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

# Additional comparisons
ATAC_fCOvmCO <- dba.report(ATAC_count,contrast = 3,th=0.05,fold=1,method = DBA_EDGER,bFlip = TRUE)
ATAC_mEBvmCO <- dba.report(ATAC_count,contrast = 4,th=0.05,fold=1,method = DBA_EDGER,bFlip = TRUE)
ATAC_fEBvmCO <- dba.report(ATAC_count,contrast = 5,th=0.05,fold=1,method = DBA_EDGER,bFlip = TRUE)
ATAC_mEBvfCO <- dba.report(ATAC_count,contrast = 6,th=0.05,fold=1,method = DBA_EDGER,bFlip = TRUE)
ATAC_fEBvfCO <- dba.report(ATAC_count,contrast = 7,th=0.05,fold=1,method = DBA_EDGER,bFlip = TRUE)
ATAC_fEBvmEB <- dba.report(ATAC_count,contrast = 8,th=0.05,fold=1,method = DBA_EDGER,bFlip = TRUE)



write.table(ATAC_fCOvmCO[ATAC_fCOvmCO$Fold>0,],"ATAC_EDGER_fCOvmCO.bed",quote = FALSE,row.names = FALSE,col.names=FALSE,sep = '\t')
write.table(ATAC_fCOvmCO[ATAC_fCOvmCO$Fold<0,],"ATAC_EDGER_mCOvfCO.bed",quote = FALSE,row.names = FALSE,col.names=FALSE,sep = '\t')
write.table(ATAC_mEBvmCO[ATAC_mEBvmCO$Fold>0,],"ATAC_EDGER_mEBvmCO.bed",quote = FALSE,row.names = FALSE,col.names=FALSE,sep = '\t')
write.table(ATAC_mEBvmCO[ATAC_mEBvmCO$Fold<0,],"ATAC_EDGER_mCOvmEB.bed",quote = FALSE,row.names = FALSE,col.names=FALSE,sep = '\t')
write.table(ATAC_fEBvmCO[ATAC_fEBvmCO$Fold>0,],"ATAC_EDGER_fEBvmCO.bed",quote = FALSE,row.names = FALSE,col.names=FALSE,sep = '\t')
write.table(ATAC_fEBvmCO[ATAC_fEBvmCO$Fold<0,],"ATAC_EDGER_mCOvfEB.bed",quote = FALSE,row.names = FALSE,col.names=FALSE,sep = '\t')
write.table(ATAC_mEBvfCO[ATAC_mEBvfCO$Fold>0,],"ATAC_EDGER_mEBvfCO.bed",quote = FALSE,row.names = FALSE,col.names=FALSE,sep = '\t')
write.table(ATAC_mEBvfCO[ATAC_mEBvfCO$Fold<0,],"ATAC_EDGER_fCOvmEB.bed",quote = FALSE,row.names = FALSE,col.names=FALSE,sep = '\t')
write.table(ATAC_fEBvfCO[ATAC_fEBvfCO$Fold>0,],"ATAC_EDGER_fEBvfCO.bed",quote = FALSE,row.names = FALSE,col.names=FALSE,sep = '\t')
write.table(ATAC_fEBvfCO[ATAC_fEBvfCO$Fold<0,],"ATAC_EDGER_fCOvfEB.bed",quote = FALSE,row.names = FALSE,col.names=FALSE,sep = '\t')
write.table(ATAC_fEBvmEB[ATAC_fEBvmEB$Fold>0,],"ATAC_EDGER_fEBvmEB.bed",quote = FALSE,row.names = FALSE,col.names=FALSE,sep = '\t')
write.table(ATAC_fEBvmEB[ATAC_fEBvmEB$Fold<0,],"ATAC_EDGER_mEBvfEB.bed",quote = FALSE,row.names = FALSE,col.names=FALSE,sep = '\t')
