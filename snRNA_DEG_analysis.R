# Load in libraries
library(plyr)
library(dplyr)
library(tibble)
library(Seurat)
library(data.table)
library(umap)
library(biomaRt)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(AnnotationFilter)
library(BiocParallel)
library(DESeq2)
library(Matrix.utils)
library(purrr)
library(GenomicFeatures)
library(tidyr)

register(BPPARAM = MulticoreParam(4))

# Read in BNST snRNA-seq 10X files
bnst_male <- Read10X(data.dir = "/BNST_snRNAseq_data/Male")
bnst_male_seurat <- CreateSeuratObject(counts = bnst_male)
bnst_female <- Read10X(data.dir = "/BNST_snRNAseq_data/Female")
bnst_female_seurat <- CreateSeuratObject(counts = bnst_female)

# Merge male and female Seurat objects
bnst_seurat <- merge(bnst_male_seurat,y = bnst_female_seurat)

# Read in metadata, add cluster, sex, and replicate as metadata
bnst_clusters_ident <- read.table(file = "BNST_metadata.csv",row.names="Barcodes",header=TRUE,sep=",")
bnst_seurat <- AddMetaData(object = bnst_seurat,metadata = bnst_clusters_ident$Cluster,col.name = 'Cluster')
bnst_seurat <- AddMetaData(object = bnst_seurat,metadata = bnst_clusters_ident$Sex,col.name = 'Sex')

sample_names_split <- as.data.frame(matrix(unlist(strsplit(rownames(bnst_clusters_ident),split = "(?=[0-9])",perl=T)),ncol=3,byrow=TRUE))
sample_names_split <- paste0(sample_names_split[,1], sample_names_split[,2])
bnst_seurat <- AddMetaData(object = bnst_seurat,metadata = sample_names_split,col.name = 'Sample')
bnst_seurat$Sample <- as.factor(bnst_seurat$Sample)
bnst_seurat$Cluster <- gsub("_", "", as.vector(bnst_seurat$Cluster)) 
bnst_seurat <- AddMetaData(object = bnst_seurat,metadata = paste0(bnst_seurat$Cluster,bnst_seurat$Sample),col.name = 'Cluster_sample')

# Transform Seurat object to SingleCellExperiment
bnst_seurat_sce <- as.SingleCellExperiment(bnst_seurat)

# Filter main Seurat object by expression and remove pseudogenes prior to DEG calling
txdb <- makeTxDbFromUCSC(genome="mm10", table="refGene")
txdb_gene_symbols <- AnnotationDbi::select(org.Mm.eg.db,
                                           keys = keys(txdb),
                                           columns = c("SYMBOL"),
                                           keytype = "ENTREZID")
bnst_seurat_sce_filter <- bnst_seurat_sce[rowSums(counts(bnst_seurat_sce) > 1) >= 5, ]
bnst_seurat_sce_filter <- bnst_seurat_sce_filter[rownames(bnst_seurat_sce_filter) %in% txdb_gene_symbols$SYMBOL,]

# Set up empty dataframes and lists to write into
df <- data.frame(male_DEGs=numeric(0),female_DEGs=numeric(0))
resOrdered_list <- list()
dds_sex_list <- list()

# Single cluster DEG analysis
for (i in unique(bnst_seurat_sce_filter$Cluster)){
  try({
  ## Drop out replicate:cluster combinations with fewer than 20 cells per cluster  
  bnst_seurat_sce_i <- bnst_seurat_sce_filter[,bnst_seurat_sce_filter$Cluster==i]
  cell_counts <- table(bnst_seurat_sce_i$Sample)
  threshold <- cell_counts[cell_counts > 20]
  bnst_seurat_sce_i_filter <- bnst_seurat_sce_i[,bnst_seurat_sce_i$Sample %in% names(threshold)]
  bnst_seurat_sce_i_filter$Sample <- droplevels(bnst_seurat_sce_i_filter$Sample)
  
  ## Aggregate counts 
  groups_i <- colData(bnst_seurat_sce_i_filter)[, c("ident", "Sample")]
  pb_i <- aggregate.Matrix(t(counts(bnst_seurat_sce_i_filter)), groupings = groups_i, fun = "sum")
  splitf_i <- sapply(stringr::str_split(rownames(pb_i), pattern = "_",  n = 2), `[`, 1)

  ## Wrangle the dataframe
  pb_i <- split.data.frame(pb_i, 
                       factor(splitf_i)) %>%
    lapply(function(u) 
      magrittr::set_colnames(t(u), 
                  stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

  ## Convert list of dataframes into single dataframe with gene names as row and sample name as column
  pb_i <- as.data.frame(do.call(cbind, pb_i))

  ## Generate metadata sheet, and ensure row order matches column order of pb object
  gg_df <- data.frame(cluster_id = i,
                      sample_id = bnst_seurat_sce_i_filter$Sample)
  gg_df$sex <- ifelse(stringr::str_detect(gg_df$sample_id,"FEMALE"),"Female", "Male")

  metadata_i <- gg_df %>%
    dplyr::select(cluster_id, sample_id, sex) 
  metadata_i <- metadata_i[!duplicated(metadata_i), ]
  metadata_i$cluster_sex <- paste0(metadata_i$cluster_id, metadata_i$sex)
  rownames(metadata_i) <- NULL
  metadata_i <- metadata_i[match(names(pb_i), metadata_i$sample_id),]

  ## Create DESeq2 object and run analysis
  dds_sex_i <- DESeqDataSetFromMatrix(pb_i,
                                    colData = metadata_i, 
                                    design = ~ cluster_sex)

  dds_sex_i <- DESeq(dds_sex_i, parallel = TRUE)
  
  ## Collect and write results
  res <- results(dds_sex_i, contrast = c("cluster_sex",paste0(i,"Male"),paste0(i,"Female")),alpha = 0.1)
  resOrdered <- res[order(res$pvalue),]
  resOrdered_na <- na.omit(resOrdered)
  write.csv(resOrdered_na, file=paste0(i,"padj.csv"))

  ## Save resOrdered to list
  resOrdered_list[[i]] <- resOrdered_na
  
  ## Save dds obj to list
  dds_sex_list[[i]] <- dds_sex_i
  
  ## Save DEG num per cluster
  gene_number <- data.frame(male_DEGs=length(rownames(resOrdered_na[(resOrdered_na$padj < 0.1) & (resOrdered_na$log2FoldChange > 0),])),female_DEGs=length(rownames(resOrdered_na[(resOrdered_na$padj < 0.1) & (resOrdered_na$log2FoldChange < 0),])))
  df[nrow(df) + 1,] = rbind(gene_number)
  rownames(df)[nrow(df)] <- paste0(i)
  })
}
write.table(df, "DEG_num_per_cluster_padj0.1.txt",sep = "\t",quote = FALSE)


# Filter for ERa+ clusters
bnst_seurat_sce_filter_Esr1 <- bnst_seurat_sce_filter[,bnst_seurat_sce_filter$Cluster %in% c("BNSTprSt18",
                                                            "BNSTpTac2",
                                                            "BNSTpBnc2",
                                                            "BNSTpEpsti1",
                                                            "BNSTprEsr2",
                                                            "BNSTpNxph2",
                                                            "BNSTpHaus4")]

# Drop out replicate-clusters pairs with fewer than 20 nuclei
cell_counts <- table(bnst_seurat_sce_filter_Esr1$Cluster_sample)
threshold <- cell_counts[cell_counts > 20]
bnst_seurat_sce_filter_Esr1_filter <- bnst_seurat_sce_filter_Esr1[,bnst_seurat_sce_filter_Esr1$Cluster_sample %in% names(threshold)]

# Aggregate counts per replicate
groups <- colData(bnst_seurat_sce_filter_Esr1_filter)[, c("ident", "Sample")]
pb <- aggregate.Matrix(t(counts(bnst_seurat_sce_filter_Esr1_filter)), groupings = groups, fun = "sum")
splitf <- sapply(stringr::str_split(rownames(pb), pattern = "_",  n = 2), `[`, 1)

# Wrangle the dataframe
pb <- split.data.frame(pb, factor(splitf)) %>%
      lapply(function(u) 
        magrittr::set_colnames(t(u), 
                               stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))
for (i in seq_along(pb)){
  colnames(pb[[i]]) <- paste0(names(pb)[i],colnames(pb[[i]]))
}


# Convert list of dataframes into single dataframe with gene names as row and sample name as column
pb <- as.data.frame(do.call(cbind, pb))

# Generate metadata sheet, and ensure row order matches column order of pb object
gg_df <- data.frame(cluster_id = bnst_seurat_sce_filter_Esr1_filter$Cluster,sample_id = bnst_seurat_sce_filter_Esr1_filter$Sample)
gg_df$sex <- ifelse(stringr::str_detect(gg_df$sample_id,"FEMALE"),"Female", "Male")
metadata <- gg_df %>% dplyr::select(cluster_id, sample_id, sex) 
metadata <- metadata[!duplicated(metadata), ]
metadata$cluster_sex <- paste0(metadata$cluster_id, metadata$sex)
metadata$cluster_sample <- paste0(metadata$cluster_id, metadata$sample_id)
rownames(metadata) <- NULL
metadata <- metadata[match(names(pb), metadata$cluster_sample),]

# Create DESeq2 object and run analysis
dds_cluster <- DESeqDataSetFromMatrix(pb,
                                    colData = metadata, 
                                    design = ~ cluster_id)

dds_cluster <- DESeq(dds_cluster,betaPrior = TRUE,parallel = TRUE)

# Collect and write results
res <- results(dds_cluster, 
               contrast = c(0,-1/6,-1/6,-1/6,-1/6,-1/6,1,-1/6),
               lfcThreshold = 2,
               altHypothesis = "greater",
               alpha = 0.01)
resOrdered <- res[order(res$pvalue),]
resOrdered_na <- na.omit(resOrdered)
write.csv(resOrdered_na, file="BNSTpr_St18_res.csv",quote = FALSE)
