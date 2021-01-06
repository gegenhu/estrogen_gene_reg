# Load in libraries
library(DESeq2)
library(ggplot2)
library(BiocParallel)
library(DiffBind)

# Load in sample sheet 
Nfix_CR <- dba(sampleSheet = "sample_sheet_peak.csv", peakFormat="narrow")
Nfix_CR_count <- dba.count(Nfix_CR,minOverlap=1,bParallel=TRUE,summits = 100)
Nfix_CR_count <- dba.contrast(Nfix_CR_count, categories=c(DBA_TREATMENT,DBA_CONDITION),minMembers = 2,block = DBA_REPLICATE)

# Run differential analysis
Nfix_CR_count <-dba.analyze(Nfix_CR_count,method=DBA_EDGER)

# Collect peaks and write to file
Nfix_CR_sex.DB <- dba.report(Nfix_CR_count,th=0.1,fold=1,method = DBA_EDGER_BLOCK,contrast = 1)
Nfix_CR_treatment.DB <- dba.report(Nfix_CR_count,th=0.1,fold=1,method = DBA_EDGER_BLOCK,contrast = 2)
write.table(Nfix_CR_sex.DB,"Nfix_EDGER_sex_peaks.bed",quote=FALSE,sep='\t')
write.table(Nfix_CR_treatment.DB,"Nfix_EDGER_treatment_peaks.bed",quote=FALSE,sep='\t')
