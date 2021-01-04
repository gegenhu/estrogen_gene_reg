# Load in library
library(DESeq2)

# Read in data
bnstData <- read.csv("adult_BNST.csv")
rownames(bnstData) <- bnstData[,1]
bnstData[,1] <- NULL

# Set up data frame
sampleNames <- c("CO_M2","CO_M3","CO_M4","CO_M5","EB_M2","EB_M3","EB_M4","EB_M5","CO_F2","CO_F3","CO_F4","CO_F5","EB_F2","EB_F3","EB_F4","EB_F5")
sampleSex <- c("male","male","male","male","male","male","male","male","female","female","female","female","female","female","female","female")
sampleBatch <- c("1","2","3","4","1","2","3","4","1","2","3","4","1","2","3","4")
sampleHormone <- c("CO","CO","CO","CO","EB","EB","EB","EB","CO","CO","CO","CO","EB","EB","EB","EB")
colData <- data.frame(sampleName=sampleNames, sex=sampleSex, hormone=sampleHormone, batch=sampleBatch)
bnstdds$hormone <- relevel(bnstdds$hormone, ref="CO")
bnstdds$sex <- relevel(bnstdds$sex, ref="female")

# one way comparison of hormone alone
bnstdds <- DESeqDataSetFromMatrix(countData = bnstData,
                                colData = colData,
                                design = ~ batch + hormone)
bnstdds <- DESeq(bnstdds)
bnst_hormone <- results(bnstdds, contrast=c("hormone","EB","CO"))
write.csv(bnst_hormone,"BNST-EBvCO.csv")

#one way comparison of sex alone
bnstdds <- DESeqDataSetFromMatrix(countData = bnstData,
                                colData = colData,
                                design = ~ batch + sex)
bnstdds <- DESeq(bnstdds)
bnst_sex <- results(bnstdds, contrast=c("sex","male","female"))
write.csv(bnst_sex,"BNST-MvF.csv")

#two way comparison of sex and hormone
bnstdds$group <- factor(paste0(bnstdds$sex, bnstdds$hormone))
design(bnstdds) <- ~ batch + group
bnstdds <- DESeq(bnstdds)

bnst_group_male <- results(bnstdds, contrast=c("group","maleEB","maleCO"))
bnst_group_female <- results(bnstdds, contrast=c("group","femaleEB","femaleCO"))
bnst_group_CO <- results(bnstdds, contrast=c("group","maleCO","femaleCO"))
bnst_group_EB <- results(bnstdds, contrast=c("group","maleEB","femaleEB"))
write.csv(bnst_group_male,"BNST-EB_MvCO_M.csv")
write.csv(bnst_group_female,"BNST-EB_FvCO_F.csv")
write.csv(bnst_group_CO,"BNST-CO_MvCO_F.csv")
write.csv(bnst_group_EB,"BNST-EB_MvEB_F.csv")

#4-way comparison of sex and hormone
design(bnstdds) <- ~ 0 + group
bnstdds <- DESeq(bnstdds)
FemCO <- results(bnstdds,contrast=c(1,-1/3,-1/3,-1/3))
FemEB <- results(bnstdds,contrast=c(-1/3,1,-1/3,-1/3))
MalCO <- results(bnstdds,contrast=c(-1/3,-1/3,1,-1/3))
MalEB <- results(bnstdds,contrast=c(-1/3,-1/3,-1/3,1))

altFemCO <- results(bnstdds,altHypothesis="greater",contrast=c(1,-1/3,-1/3,-1/3))
altFemEB <- results(bnstdds,altHypothesis="greater",contrast=c(-1/3,1,-1/3,-1/3))
altMalCO <- results(bnstdds,altHypothesis="greater",contrast=c(-1/3,-1/3,1,-1/3))
altMalEB <- results(bnstdds,altHypothesis="greater",contrast=c(-1/3,-1/3,-1/3,1))
write.csv(altFemCO,"altFemCO.csv")
write.csv(altMalEB,"altMalEB.csv")
write.csv(altFemEB,"altFemEB.csv")
write.csv(altMalCO,"altMalCO.csv")
