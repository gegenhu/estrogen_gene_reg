# Load in libraries
library(DESeq2)
library(DiffBind)
library(BiocParallel)
library(GenomicRanges)

# Load in sample sheet
ERa_CR <- dba(sampleSheet = "sample_sheet.csv", peakFormat="narrow")
ERa_CR_count <- dba.count(ERa_CR,minOverlap=1,bParallel=TRUE,summits=100)

# Establish a contrast for differential analysis
ERa_CR_count <- dba.contrast(ERa_CR_count, categories=c(DBA_CONDITION,DBA_TREATMENT),minMembers = 2,block = DBA_REPLICATE)

# Run differential analysis 
ERa_CR_count <-dba.analyze(ERa_CR_count,method=DBA_EDGER,bParallel = TRUE)

# Generate a GRanges object of the significant peaks for the contrasts
ERa_CR_count_pairwise_treatment.DB <- dba.report(ERa_CR_count,th=0.1,contrast=1,method = DBA_EDGER_BLOCK)
ERa_CR_count_pairwise_sex_EB.DB <- dba.report(ERa_CR_count,th=0.1,contrast=4,method = DBA_EDGER_BLOCK)
ERa_CR_count_pairwise_male_EB.DB <- dba.report(ERa_CR_count,th=0.1,contrast=5,method = DBA_EDGER_BLOCK)
ERa_CR_count_pairwise_female_EB.DB <- dba.report(ERa_CR_count,th=0.1,contrast=7,method = DBA_EDGER_BLOCK)

# Write data to file
write.table(ERa_CR_count_pairwise_treatment.DB,'ERa_EB_EDGER_EBvsVeh.bed',quote = FALSE,sep = '\t')
write.table(ERa_CR_count_pairwise_sex_EB.DB,'ERa_EB_EDGER_sex.bed',quote = FALSE,sep = '\t')
write.table(ERa_CR_count_pairwise_male_EB.DB,'ERa_EB_EDGER_male_EB.bed',quote = FALSE,sep = '\t')
write.table(ERa_CR_count_pairwise_female_EB.DB,'ERa_EB_EDGER_female_EB.bed',quote = FALSE,sep = '\t')
