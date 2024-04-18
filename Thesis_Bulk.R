library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)




mRNA <- read.csv("~/wyc_20231021/Bulk_SEQ/20230525测序转录组结果/mRNA.csv",row.names = 1)

mRNA <- mRNA[rowMeans(mRNA) > 1,]
group_list <- factor(c(rep("aAPCs",times=3),rep("PD1",times=3)))

Data <- data.frame(row.names = colnames(mRNA),
                   group=group_list)

dds <- DESeqDataSetFromMatrix(countData = mRNA,
                              colData = Data,
                              design = ~ group)
dds <- estimateSizeFactors(dds) 
normalized_counts <- counts(dds,normalized=T) 
dds2 <- DESeq(dds,test = "Wald")
res <- results(dds2, contrast = c("group","PD1","aAPCs"))##注意，实验组在前，对照在后，可以根据差异基因回到源表达矩阵中对照
res <- res[order(res$pvalue),]
summary(res)
my_result <- as.data.frame(res)
my_result <- na.omit(my_result)
my_result$Gene_symbol <- rownames(my_result)


####   选取p<0.2作为差异富集基因   ###
my_result <- my_result[my_result$pvalue<0.2,]

df <- bitr(unique(my_result$Gene_symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
my_result <- merge(my_result,df,by.y='SYMBOL',by.x='Gene_symbol')
data_all_sort <- my_result %>% 
  arrange(desc(log2FoldChange)) #log2FC排序,重要！！！！！！！

geneList<- data_all_sort$log2FoldChange #把foldchange按照从大到小提取出来
names(geneList) <- data_all_sort$ENTREZID #给上面提取的foldchange加上对应上ENTREZID
head(geneList)

###   kegg  ###
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 10000,
               minGSSize    = 10,
               maxGSSize    = 500,
               pvalueCutoff = 1,#
               pAdjustMethod = "none" )

kk2_symbol <- DOSE::setReadable(kk2,
                                OrgDb = org.Hs.eg.db,
                                keyType = 'ENTREZID')

kk2_symbol_result <-as.data.frame(kk2_symbol@result)
write.csv(kk2_symbol_result,file="GSEA_kegg_pathway_PD1.csv")



###  GO  ###
ego <- gseGO(geneList     =  geneList,
             OrgDb        = org.Hs.eg.db,
             ont          = "BP", pvalueCutoff = 1,  pAdjustMethod = "BH")
ego_symbol <- DOSE::setReadable(ego,
                                OrgDb = org.Hs.eg.db,
                                keyType = 'ENTREZID')
ego_symbol_result <-as.data.frame(ego_symbol@result)
write.csv(ego_symbol_result,file = "GSEA_GO.csv")



###  归一化的热图   ###
GOset <- clusterProfiler::read.gmt("~/wyc_20231021/single_cell_SEQ/define-malignant/score_methods/AUCell/c5.all.v2023.2.Hs.symbols.gmt")  #GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION
pathgene <- subset(GOset , GOset$term %in% "GOBP_CELL_CHEMOTAXIS")
pathgene <- pathgene$gene

# 筛选数据框中Gene_Symbol包含在pathgene中，且log2FoldChange>0的行
filtered_DF <- my_result %>% 
  filter(Gene_symbol %in% pathgene & log2FoldChange > 0)

# 提取符合条件的Gene_Symbol
pathgene <- filtered_DF$Gene_Symbol


pathgene <- c("HLA-DRA","HLA-DRB1","CCL5",
              "CD3D","IL2RA","STAT6","FCGR3A",
              "TNFRSF4","CCR7","HLA-DMA")

pathgene <- c("PRF1","GZMA","GZMB",
              "GNLY","GZMH","GZMK","IFNG",
              "TNF")
pathgene <- c("CXCR4","CCR5","CCR6","CXCR6")

scoregene <- readxl::read_xlsx("Score.xlsx")

pathgene <- scoregene$`Activation:Effector function`
pathgene <- na.omit(pathgene)

INT_mRNA <- mRNA[pathgene,]

new_df <- data.frame(
  ctrl_1 = rep(1, nrow(INT_mRNA)),
  ctrl_2 = rep(1, nrow(INT_mRNA)),
  ctrl_3 = rep(1, nrow(INT_mRNA)),
  aPD_1 = INT_mRNA$HD1_PD1 / INT_mRNA$HD1_aAPCs,
  aPD_2 = INT_mRNA$HD3_PD1 / INT_mRNA$HD3_aAPCs,
  aPD_3 = INT_mRNA$HD4_PD1 / INT_mRNA$HD4_aAPC,row.names = rownames(INT_mRNA)
)

heatmap_pathgene <- as.data.frame(new_df)
heatmap_pathgene <- na.omit(heatmap_pathgene)
heatmap_pathgene <- heatmap_pathgene[!apply(heatmap_pathgene, 1, function(x) any(is.infinite(x))), ]
pheatmap(heatmap_pathgene,
         scale = "row",
         color = colorRampPalette(colors = c("#177cb0","white","#ff3300"))(100),
         fontsize = 10,
         cluster_cols = F,
         cluster_rows = T,
         border=F,
         cellwidth = 18,
         cellheight = 12)























