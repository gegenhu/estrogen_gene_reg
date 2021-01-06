# Load in libraries
library(DiffBind)
library(DESeq2)
library(BiocParallel)

# Read in Franco et al. 2015 MCF-7 ChIP samples
ChIP_samples <- dba(sampleSheet='sample_sheet_ChIP.csv')
ChIP_samples_counts <- dba.count(ChIP_samples, summits=100, bParallel = T)
ChIP_samples <- dba.contrast(ChIP_samples, group1=1:2,group2=3:4,name1="E2",name2="Veh")

# Run differential analysis
ChIP_samples <- dba.analyze(ChIP_samples,method=DBA_DESEQ2,bParallel = T)

# Collect differential peaks and write to file
ChIP_results.DB <- dba.report(ChIP_samples,th=0.01,method = DBA_DESEQ2)
write.table(ChIP_results.DB[ChIP_results.DB$Fold>0,],file='MCF7_ERa_ChIP_peaks.bed',sep="\t",quote = FALSE)

# Read in MCF-7 CUT&RUN samples
CnR_samples <- dba(sampleSheet='sample_sheet_C&R.csv')
CnR_samples <- dba.count(CnR_samples, minOverlap = 1, summits=100, bParallel = TRUE)
CnR_samples <- dba.contrast(CnR_samples,categories = DBA_TREATMENT, minMembers = 2)

# Run differential analysis
CnR_samples <- dba.analyze(CnR_samples,method=DBA_DESEQ2, bParallel = TRUE)

# Collect differential peaks and write to file
CnR_results.DB <- dba.report(CnR_samples,th=0.1,method = DBA_DESEQ2)
write.table(CnR_results.DB[CnR_results.DB$Fold>0,],file='MCF7_ERa_CnR_peaks.bed',sep="\t",quote = FALSE)
