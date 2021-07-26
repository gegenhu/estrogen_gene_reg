# Load in libraries
library(DESeq2)
library(edgeR)
library(DiffBind)
library(BiocParallel)
library(GenomicRanges)

# Load in sample sheet
ERa_CR <- dba(sampleSheet = "sample_sheet_ERa_CnR.csv", peakFormat="narrow")
ERa_CR_count <- dba.count(ERa_CR,minOverlap=1,score = DBA_SCORE_READS,bParallel=TRUE,summits=100)

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

# Export peakset with raw counts for analysis with edgeR design
ERa_CR_count_df <- as.data.frame(dba.peakset(ERa_CR_count,bRetrieve=TRUE))
rownames(ERa_CR_count_df) <- paste0(ERa_CR_count_df$seqnames,'-',ERa_CR_count_df$start,'-',ERa_CR_count_df$end)
ERa_CR_count_df <- ERa_CR_count_df[c(6,7,8,9,10,11,12,13)]

# Create metadata dataframe for edgeR analysis
coldata <- data.frame(Sample = c('Male_EB_1','Male_EB_2','Male_veh_1','Male_veh_2','Female_EB_1','Female_EB_2','Female_veh_1','Female_veh_2'),
					  Treatment = c('EB','EB','Veh','Veh','EB','EB','Veh','Veh'),
					  Condition = c('Male','Male','Male','Male','Female','Female','Female','Female'),
					  Replicate = c('1','2','1','2','1','2','1','2'),
					  Group = c('Male_EB','Male_EB','Male_veh','Male_veh','Female_EB','Female_EB','Female_veh','Female_veh'))
rownames(coldata) <- coldata$Sample

# Run edgeR QLF-test with Condition:Treatment interaction design
y <- DGEList(counts=ERa_CR_count_df, samples=coldata)
y <- calcNormFactors(y)
design <- model.matrix(~Replicate + Condition + Treatment + Condition:Treatment,data=y$samples)
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=5)
res <- as.data.frame(topTags(qlf,n=Inf))
write.table(res,'ERa_CR_sex_E2_interaction.txt',quote = FALSE,sep = '\t')
