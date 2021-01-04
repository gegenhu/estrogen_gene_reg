# Load in libraries
library(DESeq2)

# Read in counts and column data
cts <- as.matrix(read.csv("Neonatal_RNA.csv",row.names="Geneid"))
coldata <- read.csv("coldata.csv", row.names=1)

# Set coldata as factors
coldata$condition <- factor(coldata$condition)
coldata$sex <- factor(coldata$sex)
coldata$sex.condition <- factor(coldata$sex.condition)

# Create DESeq object
dds <- DESeqDataSetFromMatrix(countData = cts,colData = coldata, design = ~ condition)

# Pre-filter dds object by expression
keep <- rowMeans(counts(dds)) >= 5
dds <- dds[keep,]

# Run DESeq2 and write results
dds <- DESeq(dds)
res <- results(dds,contrast = c("condition","EB","Veh"),alpha = 0.05,independentFiltering=FALSE)
resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), file="Neonatal_EBvVeh_results.csv")
