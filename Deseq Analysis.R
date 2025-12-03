library(readr)
library(dplyr)
library(DESeq2)
library(edgeR)
library(purrr)
library(readxl)
library(data.table)
library(org.Dr.eg.db)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(pheatmap)
library(ggplot2)
library(enrichplot)


count_data <- read_table("/Users/ardicsahin/Desktop/Bioinformatics/count_data_21days.tabular")
count_data <- as.data.frame(count_data)

count_data <- count_data[-c(1),-c(8:55) ]
colnames(count_data) <- c("GENEID","CMB_2", "CMB_1", "NORM_2", "NORM_1", "TrA_2", "TrA_1")
count_data$TrA_2 <- as.character(count_data$TrA_2)
count_data <- count_data[!duplicated(count_data$GENEID), ]

GeneID_col <- count_data[, -c(2:7)]
GeneID_col <- as.data.frame(GeneID_col)
colnames(GeneID_col) <- c("GENEID")

rownames(count_data) <-count_data$GENEID
count_data <- count_data[, -1]

count_data <- na.omit(count_data)
count_data <- data.frame(lapply(count_data, as.numeric))

count_dataCvN <- count_data[,-c(1,2)]
count_dataTvC <- count_data[, -c(3,4)]
count_dataTvN <- count_data[, -c(5,6)]

#All Groups

#TrAvsNORM

conditions <- c(rep("CMB",2),rep("TrA",2),rep("NORM",2))

sample.design <- data.frame(condition= conditions, row.names = colnames(count_data))
count_data <- na.omit(count_data)

ddset <- DESeqDataSetFromMatrix(countData = count_data, colData = sample.design, design =~ condition)
dds <- estimateSizeFactors(ddset)
dds <- ddset[rowSums(counts(dds))>0,]
dds <- DESeq(dds)

TrAvsNORM <- data.frame(results(dds, contrast= c("condition", "TrA", "NORM")))
CMBvsTrA <- data.frame(results(dds, contrast= c("condition", "CMB", "TrA")))
CMBvsNORM <- data.frame(results(dds, contrast= c("condition", "CMB", "NORM")))

TrAvsNORM$GENEID <- row.names(TrAvsNORM)
CMBvsTrA$GENEID <- row.names(CMBvsTrA)
CMBvsNORM$GENEID <- row.names(CMBvsNORM)

names(TrAvsNORM)[-7] <- paste0("TrAvsNORM_", names(TrAvsNORM)[-7])
names(CMBvsTrA)[-7] <- paste0("CMBvsTrA_", names(CMBvsTrA)[-7])
names(CMBvsNORM)[-7] <- paste0("CMBvsNORM_", names(CMBvsNORM)[-7])

CvN_TvCmerge <- merge(TrAvsNORM,CMBvsTrA, by.x= "GENEID", by.y="GENEID", all.x=FALSE, all.y=FALSE)
CvN_TvC_TvN_merge <- merge(CvN_TvCmerge,CMBvsNORM, by.x= "GENEID", by.y="GENEID", all.x=FALSE, all.y=FALSE)
CvN_TvC_TvN_merge <- na.omit(CvN_TvC_TvN_merge)

TrA_TvC_TvN_filter <- CvN_TvC_TvN_merge[, c("GENEID","TrAvsNORM_log2FoldChange","TrAvsNORM_padj",
                                       "CMBvsTrA_log2FoldChange","CMBvsTrA_padj",
                                       "CMBvsNORM_log2FoldChange","CMBvsNORM_padj")]

CMBvsTrAfilter <- TrA_TvC_TvN_filter %>% filter(TrA_TvC_TvN_filter$CMBvsTrA_padj < 0.05) #Sadece CMBvsTrA padj significant genlerin 3 gruptaki LogFC değerlerini almak için
TrAvsNORMfilter <- TrA_TvC_TvN_filter %>% filter(TrA_TvC_TvN_filter$TrAvsNORM_padj < 0.05)
CMBvsNORMfilter <- TrA_TvC_TvN_filter %>% filter(TrA_TvC_TvN_filter$CMBvsNORM_padj < 0.05)
CMBvsTrA_LogFC <- CMBvsTrAfilter[,c("TrAvsNORM_log2FoldChange","CMBvsTrA_log2FoldChange", "CMBvsNORM_log2FoldChange")]
TrAvsNORM_LogFC <- TrAvsNORMfilter[,c("TrAvsNORM_log2FoldChange","CMBvsTrA_log2FoldChange", "CMBvsNORM_log2FoldChange")]
CMBvsNORM_LogFC <- CMBvsNORMfilter[,c("TrAvsNORM_log2FoldChange","CMBvsTrA_log2FoldChange", "CMBvsNORM_log2FoldChange")]


resTvC <- results(dds, contrast= c("condition", "CMB", "TrA")) 
resCvN <- results(dds, contrast= c("condition", "TrA", "NORM")) 
resTvN <- results(dds, contrast= c("condition", "CMB", "NORM")) 

resultsNames(dds)
pheatmap(CMBvsTrA_LogFC, scale= "row", show_rownames= T, cluster_cols= T, main="Heatmap of Significant Genes for CMB vs TrA", angle_col = c("0"))
pheatmap(TrAvsNORM_LogFC, scale= "row", show_rownames= T, cluster_cols= T, main="Heatmap of Significant Genes for TrA vs NORM", angle_col = c("0"))
pheatmap(CMBvsNORM_LogFC, scale= "row", show_rownames= T, cluster_cols= T, main="Heatmap of Significant Genes for CMB vs NORM", angle_col = c("0"))




res_ordTvC <- resTvC[order(resTvC$padj),]
top_genesTvC <- row.names(res_ordTvC)[1:20]
counts_tTvC <- counts(dds, normalized= T)
counts_topTvC <- counts_tTvC[top_genesTvC,]
log_counts_topTvC <- log2(counts_topTvC+1)
pheatmap(log_counts_topTvC, scale= "row", main="Top 20 Genes for CMBvsTrA", angle_col = c("0"))

res_ordCvN <- resCvN[order(resCvN$padj),]
top_genesCvN <- row.names(res_ordCvN)[1:20]
counts_tCvN <- counts(dds, normalized= T)
counts_topCvN <- counts_tCvN[top_genesCvN,]
log_counts_topCvN <- log2(counts_topCvN+1)
pheatmap(log_counts_topCvN, scale= "row", main="Top 20 Genes for TrAvsNORM", angle_col = c("0"))

res_ordTvN <- resTvN[order(resTvN$padj),]
top_genesTvN <- row.names(res_ordTvN)[1:20]
counts_tTvN <- counts(dds, normalized= T)
counts_topTvN <- counts_tTvN[top_genesTvN,]
log_counts_topTvN <- log2(counts_topTvN+1)
pheatmap(log_counts_topTvN, scale= "row", main="Top 20 Genes for CMBvsNORM", angle_col = c("0"))



#Enrichment


ogl_TvC <- CMBvsTrAfilter$CMBvsTrA_log2FoldChange #LogFC of CMBvsTrA
names(ogl_TvC) <- CMBvsTrAfilter$GENEID

gl_TvC <- na.omit(ogl_TvC) 


sig_TvC = subset(CMBvsTrAfilter, CMBvsTrA_padj< 0.05)
rownames(sig_TvC) <- sig_TvC$GENEID
genes_TvC <- sig_TvC$GENEID

go_enrichTvC <- enrichGO(gene = genes_TvC, OrgDb = org.Dr.eg.db,keyType = "ENTREZID", readable = TRUE, 
                        ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.10)

dotplot(go_enrichTvC)
upsetplot(go_enrichTvC)
barplot(go_enrichTvC, drop = TRUE, showCategory = 8, title = "CMBvsTrA Enrichment Data", font.size = 10)


ogl_CvN <- TrAvsNORMfilter$TrAvsNORM_log2FoldChange #LogFC of TrAvsNORM
names(ogl_CvN) <- TrAvsNORMfilter$GENEID
gl_CvN <- na.omit(ogl_CvN)

sig_CvN = subset(TrAvsNORMfilter, TrAvsNORM_padj< 0.05)
rownames(sig_CvN) <- sig_CvN$GENEID
genes_CvN <- sig_CvN$GENEID

go_enrichCvN <- enrichGO(gene = genes_CvN, OrgDb = org.Dr.eg.db,keyType = "ENTREZID", readable = TRUE, 
                        ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.10)

dotplot(go_enrichCvN)
upsetplot(go_enrichCvN)
barplot(go_enrichCvN, drop = TRUE, showCategory = 8, title = "TrAvsNORM Enrichment Data", font.size = 10)


ogl_TvN <- CMBvsNORMfilter$CMBvsNORM_log2FoldChange 
names(ogl_TvN) <- CMBvsNORMfilter$GENEID
gl_TvN <- na.omit(ogl_TvN)

sig_TvN = subset(CMBvsNORMfilter, CMBvsNORM_padj< 0.05)
rownames(sig_TvN) <- sig_TvN$GENEID
genes_TvN <- sig_TvN$GENEID

go_enrichTvN <- enrichGO(gene = genes_TvN, OrgDb = org.Dr.eg.db,keyType = "ENTREZID", readable = TRUE, 
                        ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.10)


dotplot(go_enrichTvN)
upsetplot(go_enrichTvN)
barplot(go_enrichTvN, drop = TRUE, showCategory = 8, title = "CMBvsNORM Enrichment Data", font.size = 10)
























  
  
  
  
