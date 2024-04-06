################################################################################
library(DESeq2)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(limma)
library(dplyr)
library(ggrepel)
library(ggsci)
library(clusterProfiler)
library(enrichplot)
lapply(c('clusterProfiler','enrichplot','patchwork'), function(x) {library(x, character.only = T)})
library(org.Hs.eg.db)
library(patchwork)
library(WGCNA)
library(GSEABase)
library(GSVA)
library(scatterplot3d)
library(factoextra)
library(stringr)
library(stringi)
library(biomaRt)
library(pals)
library(RColorBrewer)
library(cowplot)
library(xCell)
library(GenomicFeatures)#count转fpkm要用
# 读取数据+差异表达分析

data <- read.table("NC_Bulk_Counts_Final.txt",header = T)
data <- data[, !grepl("gene_ENSG\\.", colnames(data))]
data$gene_ENSG <- stri_sub(data$gene_ENSG,1,15)
ENSGID <- data$gene_ENSG
gene_symbol <- bitr(ENSGID, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")#可同时转两种，symbol and enterz
data <- data.frame(gene_symbol,data[match(gene_symbol$ENSEMBL,data$gene_ENSG),])
data <- data[,-c(1,3)]

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
gene_biotype <- getBM(filters = "hgnc_symbol", 
                         attributes = c("hgnc_symbol", "gene_biotype"),
                         values = data$SYMBOL, mart = mart)
colnames(gene_biotype)[which(colnames(gene_biotype) == "hgnc_symbol")] <- "SYMBOL"

data <- merge(data,gene_biotype,"SYMBOL")

anyDuplicated(data)#如果返回数字>0，说明有重复
data <- data[!duplicated(data), ]

mRNA <- data[data$gene_biotype=="protein_coding",]
mRNA <- mRNA[,-c(28)]
mRNA$mean <- rowMeans(mRNA[,2:27])
mRNA <- arrange(mRNA,desc(mean))#按mean降序排列
mRNA <- mRNA %>% distinct(SYMBOL, .keep_all = TRUE)#删除重复基因
mRNA <- dplyr::select(mRNA,-c(mean))
rownames(mRNA) <- mRNA$SYMBOL
mRNA <- mRNA[,-1]

groupall <- readxl::read_xlsx("group.xlsx")
groupall <- groupall$...1

colnames(mRNA) <- groupall

#count转fpkm数据----------------------------------------------------------------
txdb <- makeTxDbFromGFF("gencode.v36.annotation.gtf",format = 'gtf')
exons_gene <- exonsBy(txdb,by = 'gene')
exons_gene_lens <- lapply(exons_gene,function(x){sum(width(reduce(x)))})

class(exons_gene_lens)
length(exons_gene_lens)

exons_gene_lens1 <- as.data.frame(exons_gene_lens)
exons_gene_lens1 <- t(exons_gene_lens1)
exons_gene_lens1 <- as.data.frame(exons_gene_lens1)

counts <- read.table("NC_Bulk_Counts_Final.txt",header = T)
counts <- counts[, !grepl("gene_ENSG\\.", colnames(counts))]
counts$gene_ENSG <- stri_sub(counts$gene_ENSG,1,15)

xc <- gsub("\\.(\\.?\\d*)","",rownames(exons_gene_lens1))
xc
rownames(exons_gene_lens1) = xc
colnames(exons_gene_lens1) = "Length"
exons_gene_lens1$ENSG <- rownames(exons_gene_lens1)

count_with_length <- merge(counts,exons_gene_lens1,by.x= "gene_ENSG",by.y= "ENSG")
count_with_length <- count_with_length %>% distinct(gene_ENSG, .keep_all = TRUE)#删除重复基因
rownames(count_with_length) <- count_with_length$gene_ENSG
count_with_length <- count_with_length[,-1]

#先把length除以1000，单位要kb
kb <- count_with_length$Length/1000

countdata <- count_with_length[,1:26]
rpk <- countdata/kb

#FPKM计算
fpkm <- t(t(rpk)/colSums(countdata) * 10^6)
head(fpkm,5)
fpkm <- as.data.frame(fpkm)
fpkm$gene_ENSG <- rownames(fpkm)
ENSGID <- fpkm$gene_ENSG
gene_symbol <- bitr(ENSGID, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")#可同时转两种，symbol and enterz

fpkm <- merge(fpkm,gene_symbol,by.x = "gene_ENSG",by.y = "ENSEMBL")
fpkm <- fpkm[,-c(1)]

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))#可能需要尝试多遍，服务器不稳定
gene_biotype <- getBM(filters = "hgnc_symbol", 
                      attributes = c("hgnc_symbol", "gene_biotype"),
                      values = fpkm$SYMBOL, mart = mart)
colnames(gene_biotype)[which(colnames(gene_biotype) == "hgnc_symbol")] <- "SYMBOL"

fpkm <- merge(fpkm,gene_biotype,"SYMBOL")


fpkm <- fpkm %>% distinct(SYMBOL, .keep_all = TRUE)#删除重复基因
rownames(fpkm) <- fpkm$SYMBOL

anyDuplicated(fpkm)#如果返回数字>0，说明有重复
fpkm <- fpkm[!duplicated(fpkm), ]

mRNA <- fpkm[fpkm$gene_biotype=="protein_coding",]
mRNA <- mRNA[,-c(28)]
mRNA$mean <- rowMeans(mRNA[,2:27])
mRNA <- arrange(mRNA,desc(mean))#按mean降序排列
mRNA <- mRNA %>% distinct(SYMBOL, .keep_all = TRUE)#删除重复基因
mRNA <- dplyr::select(mRNA,-c(mean))
rownames(mRNA) <- mRNA$SYMBOL
mRNA <- mRNA[,-1]

groupall <- readxl::read_xlsx("group.xlsx")
groupall <- groupall$...1

colnames(mRNA) <- groupall
#保存FPKM矩阵
write.csv(mRNA,file="mRNA_fpkm.csv",sep="\t",quote=F)

################################################################################


data <- read.csv("PD1_aAPCs.csv", header = T, row.names = 1)

expr_count <- cbind(gene_type=data$gene_biotype,gene_name=data$gene_name,data)#整合
expr_count <- expr_count[,0:8]
mRNA <- expr_count[expr_count$gene_type=="protein_coding",]
mRNA <- mRNA[,-1]
dim(mRNA)
mRNA$mean <- rowMeans(mRNA[,2:6])#计算每行表达平均值
mRNA <- arrange(mRNA,desc(mean))#按mean降序排列
mRNA <- mRNA %>% distinct(gene_name, .keep_all = TRUE)#删除重复基因
mRNA <- dplyr::select(mRNA,-c(mean))

rownames(mRNA) <- mRNA$gene_name
mRNA <- mRNA[,-1]

write.csv(mRNA,file = "mRNA.csv")

#读入mRNA数据-------------------------------------------------------------------
mRNA <- read.csv("mRNA.csv",header = T)
rownames(mRNA) <- mRNA$X
mRNA <- mRNA[,-c(1)]

expr <- mRNA
expr <- expr[rowMeans(expr) > 1,]#这个1”是想要的阈值，可以筛

#创建分组因子变量的几种方法
group_list <- factor(c(rep("aAPCs",times=2),rep("PD1",times=2)))#创建分组因子变量

group_list <- ifelse(grepl("pre$", colnames(expr)),
                     "pre", "post")

group_list <- ifelse(as.numeric(str_sub(rownames(data),14,15))<10,'tumor','normal')

group_list <- factor(MPR)



Data <- data.frame(row.names = colnames(expr),
                   group=group_list)

dds <- DESeqDataSetFromMatrix(countData = expr,
                              colData = Data,
                              design = ~ group)
dds2 <- DESeq(dds)
res <- results(dds2, contrast = c("group","MPR","nonMPR"))##注意，实验组在前，对照在后，可以根据差异基因回到源表达矩阵中对照
res <- res[order(res$pvalue),]
summary(res)
my_result <- as.data.frame(res)
my_result <- na.omit(my_result)
my_result$Gene_symbol <- rownames(my_result)
my_result <- my_result %>% dplyr::select('Gene_symbol', 
                                         colnames(my_result)[1:dim(my_result)[2]-1],
                                         everything())
rownames(my_result) <- NULL

my_result$regulate <- ifelse(my_result$pvalue > 0.05, "unchanged",
                             ifelse(my_result$log2FoldChange > 1.5, "up-regulated",
                                    ifelse(my_result$log2FoldChange < -1.5, "down-regulated","unchanged")))
table(my_result$regulate)
#将上下调基因单独存储
#配对的差异基因-----------------------------------------------------------------
## 4.1表达矩阵
###先得得到mRNA###
mRNA <- mRNA[,c(1,2,4,5)]
mRNA <- mRNA[rowMeans(mRNA)>0,] # 剔除表达量低的基因
data <- mRNA
group_list <- c(rep("aAPCs",2), rep("PD-1",2)) 
group_list <- factor(group_list,levels = c("aAPCs","PD-1"))
data = apply(data, 2, as.integer) ## DESeq2分析需要是整数
data <- as.data.frame(data)
row.names(data) <- row.names(mRNA)

## 4.2分组矩阵，配对分析与常规分析最大的区别就在分组矩阵

condition = group_list
# 配对分析要加上这段代码，知道谁和谁是一对，比如1,1是一对，5,5是一对
subject <- factor(c(1,1,2,2))  

coldata <- data.frame(row.names = colnames(data), condition)

# 注意在design中加上配对信息
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = coldata,
                              design = ~subject +condition) 

dds$condition<- relevel(dds$condition, ref = "control") 

## 4.3差异表达矩阵，还是和常规分析一样

dds <- DESeq(dds)
nrDEG_DESeq2 <- as.data.frame(results(dds))
rld <- rlog(dds)
# 这里我还提取了标准化后的表达矩阵，可以用于后续的热图绘制等等
normal_gset <- assay(rld) 
nrDEG_DESeq2 = nrDEG_DESeq2[order(nrDEG_DESeq2$log2FoldChange),] 

## 4.4定义差异基因

nrDEG <- nrDEG_DESeq2
nrDEG$Group = "notsignificant"

logFC_cutoff <- 0.6
nrDEG$Group[which( (nrDEG$pvalue < 0.05) & (nrDEG$log2FoldChange > logFC_cutoff) )] = "upregulated"
nrDEG$Group[which( (nrDEG$pvalue < 0.05) & (nrDEG$log2FoldChange < -logFC_cutoff) )] = "downregulated"

table(nrDEG$Group)


#清空环境以便后续操作-----------------------------------------------------------
write.csv(my_result,file = "my_result_MPR.csv")

ls()
rm(list=ls())
gc()

my_result <- read.csv("my_result")
my_result <- my_result[,-1]
my_result <- my_result %>%
  arrange(desc(log2FoldChange))


#火山图绘制---------------------------------------------------------------------
DEG_deseq2 <- subset(my_result, pvalue < 0.2 & abs(log2FoldChange) > 0.5)#######注意要改
upgene <- DEG_deseq2[DEG_deseq2$regulate=='up-regulated',]
downgene <- DEG_deseq2[DEG_deseq2$regulate=='down-regulated',]
write.csv(DEG_deseq2,file = "DEG_PD1_AAPCS.csv")

plot(my_result$log2FoldChange,-log2(my_result$pvalue))

P <- ggplot(data=my_result, aes(x = log2FoldChange, y = -log10(pvalue), color=regulate)) +
  geom_point(shape=16, size=3) +
  theme_set(theme_set(theme_bw(base_size = 20))) +
  xlab("log2 FoldChange") +
  ylab("-log10 p-value") +
  theme(plot.title = element_text(size = 10, hjust = 2.5)) +
  theme_classic() +
  scale_color_manual(values = c('#86D47D','#DEF0DF','#C43B99')) +
  geom_vline(xintercept = c(-1.5,1.5), lty=4, col='gray', lwd=0.8) +#logFC分界???
  geom_hline(yintercept = -log10(0.05), lty=2, col='gray', lwd=0.6) +#adj.p.val分界???
  labs(title = 'MPR-nonMPR') +
  geom_text_repel(data = filter(my_result, abs(my_result$log2FoldChange) > 1.5 & -log10(my_result$pvalue) > 2),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 10),
                  aes(label = Gene_symbol, ),min.segment.length = 0,force = 25,
                  size = 5) 
plot(P)

annotate("text",x = upgene$log2FoldChange[1:10],y = (-log10(upgene$pvalue[1:10])),label=upgene$Gene_symbol[1:10],size=3.5) +
  annotate("text",x = downgene$log2FoldChange[1:10],y = (-log10(downgene$pvalue[1:10])),label=downgene$Gene_symbol[1:10],size=3.5)

#取感兴趣的通路的基因做火山图  需要用到scoregene里面的基因 my_results数据-----
scoregene <- readxl::read_xlsx("score.xlsx")
maker <- rownames(heatmap_pathgene)

rownames(my_result) <- my_result$Gene_symbol
###可以改一下regulate的cutoff条件
my_result$regulate <- ifelse(my_result$pvalue > 0.2, "unchanged",
                             ifelse(my_result$log2FoldChange > 0.2, "up-regulated",
                                    ifelse(my_result$log2FoldChange < -0.2, "down-regulated","unchanged")))

PWYL <- scoregene$Cytotoxicity
PWYL <- na.omit(PWYL)
PWYL_gene <- my_result[maker,]
PWYL_gene <- na.omit(PWYL_gene)
PWYL_gene$Gene_symbol <- row.names(PWYL_gene)

#title！show指定gene/show界限gene
P <- ggplot(data=PWYL_gene, aes(x = log2FoldChange, y = -log10(pvalue), color=regulate)) +
  geom_point(shape=16, size=3) +
  theme_set(theme_set(theme_bw(base_size = 20))) +
  xlab("log2 FoldChange") +
  ylab("-log10 p-value") +
  theme(plot.title = element_text(size = 10, hjust = 2.5)) +
  theme_classic() +
  scale_color_manual(values = c('#177cb0','gray','#ff3300')) +  #'#177cb0','gray','#ff3300'
  geom_vline(xintercept = c(-0.2,0.2), lty=4, col='gray', lwd=0.8) +#logFC分界???
  geom_hline(yintercept = -log10(0.2), lty=2, col='gray', lwd=0.6) +#adj.p.val分界???
  labs(title = 'function') +
  geom_text_repel(data = filter(PWYL_gene, Gene_symbol %in% maker),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 10),
                  aes(label = Gene_symbol, ),min.segment.length = 0,force = 10,
                  size = 5)+
  coord_cartesian(xlim=c(-max(abs(PWYL_gene$log2FoldChange)), max(abs(PWYL_gene$log2FoldChange))))
plot(P)

#show-指定的GENE:  Gene_symbol %in% c("GZMA","GZMB","GZMK","GZMH","GNLY")
#show-指定GENE界限：abs(PWYL_gene$log2FoldChange) > 0.2 & -log10(PWYL_gene$pvalue) > 0.5 
#GO富集分析---------------------------------------------------------------------


diff_enrich = my_result[(my_result$pvalue< 0.05 & abs(my_result$log2FoldChange)>0.5),]                    
diff_enrich = diff_enrich[order(diff_enrich$log2FoldChange),]


UP_enrich  <- diff_enrich[diff_enrich$log2FoldChange>0, ]
rownames(UP_enrich) <- UP_enrich[,1]

UP_enrich <- bitr(rownames(UP_enrich), fromType="SYMBOL", 
                  toType="ENTREZID",
                  OrgDb="org.Hs.eg.db")

ego_up <- enrichGO(gene = UP_enrich$ENTREZID,
                   OrgDb = org.Hs.eg.db, 
                   ont="all",
                   readable=TRUE,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.2)#只选择BP生物过程富集,CC细胞组成，MF分子功能

frame_ego_up <- as.data.frame(ego_up)
write.csv(frame_ego_up,file="GO_UP_PD1.csv")
#########-----NES值看整条通路，通路---柱状图，GSEA通路整体基因差异，KEGG看分子通路（增殖活化功能stat3）基因热图（功能活化增值趋化）
#########-----PCA分群


DOWN_enrich  <- diff_enrich[diff_enrich$log2FoldChange<0, ]
rownames(DOWN_enrich) <- DOWN_enrich[,1]

DOWN_enrich <- bitr(rownames(DOWN_enrich), fromType="SYMBOL", 
                  toType="ENTREZID",
                  OrgDb="org.Hs.eg.db")

ego_down <- enrichGO(gene = DOWN_enrich$ENTREZID,
                   OrgDb = org.Hs.eg.db, 
                   ont="all",
                   readable=TRUE,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.2)#只选择BP生物过程富集,CC细胞组成，MF分子功能

#绘制柱状???
pdf(file = "PD1-aAPCs-UP_柱状图.pdf", width = 10, height = 9)
barplot(ego_up,
        drop= TRUE,
        showCategory = 10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~.,scales = 'free')
dev.off()

#绘制气泡???
pdf(file = "PD1-aPACs-UP_气泡图.pdf", width = 10, height = 9)
dotplot(ego_up,showCategory = 10)
dev.off()

#网格图（以BP为例???
BP_EGO_UP <- enrichGO(gene = UP_enrich$ENTREZID,
                      keyType = "ENTREZID",
                      OrgDb = org.Hs.eg.db, 
                      ont="BP",
                      readable=TRUE,
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.1)
pdf("PD1-aPACs-UP_网格图.pdf",width = 10,height = 20)
plotGOgraph(BP_EGO_UP)
dev.off()

#过滤简化函数：BP_EGO_UP <- simplify(ego_up, cutoff=0.05, by="p.adjust", select_fun = min)

#网络???
BP_EGO_UP <- pairwise_termsim(BP_EGO_UP)
pdf("PD1-aAPCs-UP-BP_网络图.pdf",width = 9, height = 9)
emapplot(BP_EGO_UP,cex.params = list(category_label = .4,line = .2)) +  
  scale_color_continuous(low = "#e06663", high = "#327eba", 
                         name = "p.adjust", 
                         guide = guide_colorbar(reverse = TRUE, order=1), 
                         trans='log10')
dev.off()

#GSEA富集分析-------------------------------------------------------------------

deg <- read.csv("my_result")
deg <- deg[,-1]
my_result <- read.csv("DEG_result_epi_time.csv")
rownames(my_result) <- my_result$Gene_symbol



df <- bitr(unique(my_result$Gene_symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
DEG <- my_result
DEG <- merge(DEG,df,by.y='SYMBOL',by.x='Gene_symbol')
data_all_sort <- DEG %>% 
  arrange(desc(log2FoldChange)) #log2FC排序,重要！！！！！！！

geneList<- data_all_sort$log2FoldChange #把foldchange按照从大到小提取出来
names(geneList) <- data_all_sort$ENTREZID #给上面提取的foldchange加上对应上ENTREZID
head(geneList)

#开始kegg富集分析
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
kk2_symbol_result <-as.data.frame(kk2_symbol@result)#将kk2里面的ENTREZID转为symbol
#将感兴趣的通路的基因列为向量

CKCKR <- kk2_symbol@geneSets[["hsa04060"]]
CKCKR <- bitr(CKCKR,fromType = "ENTREZID",
              toType = "SYMBOL",
              OrgDb = org.Hs.eg.db)
NKR <- kk2_symbol@geneSets[["hsa04650"]]
NKR <- bitr(NKR,fromType = "ENTREZID",
              toType = "SYMBOL",
              OrgDb = org.Hs.eg.db)
class(kk2)
colnames(kk2@result)
kegg_result <- as.data.frame(kk2)
rownames(kk2@result)[head(order(kk2@result$enrichmentScore))]
af <- as.data.frame(kk2@result)
write.csv(kk2_symbol_result,file="GSEA_kegg_pathway_PD1.csv")

#GSEA GO------------------------------------------------------------------------
ego <- gseGO(geneList     =  geneList,
             OrgDb        = org.Hs.eg.db,
             ont          = "BP", pvalueCutoff = 1,  pAdjustMethod = "BH")
ego_symbol <- DOSE::setReadable(ego,
                                OrgDb = org.Hs.eg.db,
                                keyType = 'ENTREZID')
ego_symbol_result <-as.data.frame(ego_symbol@result)
frame_ego <- as.data.frame(ego@result)
write.csv(ego_symbol_result,file = "GSEA_GO_pathway_epi_RECIST.csv")

#感兴趣的通路单独看热图
mRNA <- mRNA_count

leukocyte_proliferation<- ego_symbol@geneSets[["GO:0070661"]]
leukocyte_proliferation <- bitr(leukocyte_proliferation,fromType = "ENTREZID",
              toType = "SYMBOL",
              OrgDb = org.Hs.eg.db)
leukocyte_gene <- mRNA[leukocyte_proliferation$SYMBOL,]
leukocyte_gene <- leukocyte_gene[which(rowSums(leukocyte_gene)>0),]#校准时一行的值不能全相同
pheatmap(leukocyte_gene,
         scale = "row",
         cluster_cols = F,
         cluster_rows = T,
         angle_col = 45,border=F,show_rownames = F,
         cutree_cols = 2,
         cutree_rows = 1)

#通路柱状图---------------------------------------------------------------------
rownames(ego_symbol_result) <- ego_symbol_result$ID
my_pathway <- c("GO:0042098","GO:0042129","GO:0090068","GO:0022402","GO:0012501",
                "GO:0008219","GO:0060326","GO:1990869","GO:0072678")

my_pathway <- c("GO:0042098","GO:0042110",
                "GO:0008219","GO:0060326","GO:1990869","GO:0072678")

my_pathway <- c("GO:0070661","GO:0046651","GO:0042102")#proliferation related

my_pathway <- ego_symbol_result[my_pathway,]
ggplot(my_pathway,aes(x=my_pathway$NES,
                      y=reorder(my_pathway$Description,my_pathway$NES)))+
  theme_bw()+
  geom_col(aes(fill=my_pathway$pvalue))+
  labs(x= "NES", y= "GO pathways")


#kegg pathway-柱状图
mykegg_pathway <- af[c("hsa04110","hsa04670","hsa04060","hsa04657","hsa04064",
                   "hsa04668","hsa04014","hsa04151","hsa04010","hsa04630",
                   "hsa05235","hsa04650","hsa04350"),] 
mykegg_pathway <- af[c("hsa04668","hsa04650"),] #function,CK
mykegg_pathway <- mykegg_pathway[,-1]

ggplot(mykegg_pathway,aes(x=mykegg_pathway$NES,
                 y=reorder(mykegg_pathway$Description,mykegg_pathway$NES)))+
  geom_col(aes(fill=mykegg_pathway$pvalue))+
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"))+
  labs(x= "NES", y= "KEGG pathways")

#kegg pathway heatmap
keggckckr <- kk2_symbol@geneSets[["hsa04350"]]
keggckckr <- bitr(keggckckr,fromType = "ENTREZID",
                                toType = "SYMBOL",
                                OrgDb = org.Hs.eg.db)
keggckckr_gene <- mRNA[keggckckr$SYMBOL,]
keggckckr_gene <- keggckckr_gene[which(rowSums(keggckckr_gene)>0),]#校准时一行的值不能全相同
pheatmap(keggckckr_gene,
         scale = "row",
         cluster_cols = F,
         cluster_rows = T,
         angle_col = 45,border=F,show_rownames = F,
         cutree_cols = 2,
         cutree_rows = 1)


#排序后取前5个和后5个一起展示
num=5
pdf(paste0("2.","all_GSEA.pdf"),width = 10,height = 10)
gseaplot2(kk2, geneSetID = rownames(kk2@result)[c(head(order(kk2@result$enrichmentScore),num),tail(order(kk2@result$enrichmentScore),num))])
dev.off()

#GSEA-PLOT----------------------------------------------------------------------
enrichplot::gseaplot2(kk2,
          title = "pathway",  #设置标题
          c('hsa04630','hsa04151'), #绘制hsa04658通路的结果，通路名称与编号对应 #线条颜色
          base_size = 6, #基础字体的大小
          subplots = 1:3, 
          rel_heights = c(1.5, 0.5, 1),
          pvalue_table = T) # 显示p值????????????


enrichplot::gseaplot2(kk2,
                      c('hsa04151','hsa04010','hsa04630'),
                      base_size = 17,
                      pvalue_table = T,title = 'KEGG-GSEA')

#山脊图，展示10个，自行保存
kk2@result$Description=gsub("HALLMARK_","",kk2@result$Description)
ridgeplot(kk2,showCategory = 20)

kegg <- enrichKEGG(df$ENTREZID,
                   organism = "hsa",
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)
write.table(kegg$ID, file = "KEGG_IDs.txt", #将所有KEGG富集到的通路写入本地文件查看
            append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
browseKEGG(kegg,"hsa04060")


#kegg-pathview------------------------------------------------------------------
#提取一个含ENTERZID和log2fc的矩阵
#首先需要DEG中的ENTERZID&LOG2FC
library(pathview)
library(tidyverse)
kgp <- DEG
row.names(kgp) <- kgp$ENTREZID
kgp2 <- dplyr::select(kgp,log2FoldChange)
pathview(gene.data = kgp2, #上面包含行名为entrezID的logFC值的数据框
         pathway.id = "hsa04060", #选择一个KEGG信号通路
         species = "hsa",
         out.suffix = "3") #图的后缀
#免疫浸润xcell分析--------------------------------------------------------------
genes <- xCell.data$genes
signatures <- xCell.data$signatures
signatures <- as.data.frame(signatures)
#计算富集分数
scores = rawEnrichmentAnalysis(mRNA, signatures,genes, 
                                parallel.sz=4, parallel.type = "SOCK")
#富集分数转换
tscores <-  transformScores(scores,
                            fit.vals = xCell.data$spill$fv,
                            scale = T)
#溢出补偿
result <- spillOver(tscores,K = xCell.data$spill$K, alpha = 0.5)

#一步法
result <- xCellAnalysis(mRNA)

#可视化处理
result <- as.data.frame(result)

#免疫浸润IOBR分析（整合）
library(IOBR)
library(tidyverse)
library(reshape2)
#查看计算方法
tme_deconvolution_methods
#返回特征估计所用的参数
signature_score_calculation_methods
#相关基因集--tme  metabolism  metabolism  collection
signature_tme
signature_metabolism
signature_metabolism
signature_collection

#cibersort#
cibersort <- deconvo_tme(eset= mRNA,
                         method = "cibersort",
                         arrays = F,
                         perm = 10#最好大于100次
                         )

cell_bar_plot(input = cibersort,
                    title = "Cell Fractin",
                    legend.position = "bottom",
                    palette = 3,
                    show_col = F,
                    coord_filp = T)

xcibersort <- cibersort
cibersort <- cibersort[c(1,2,7,13,14,20,24),]

#转换为长数据，方便画图
cibersort_long <- cibersort %>% 
  select(`P-value_CIBERSORT`,Correlation_CIBERSORT, RMSE_CIBERSORT,ID,everything()) %>% 
  pivot_longer(- c(1:4),names_to = "cell_type",values_to = "fraction") %>% 
  dplyr::mutate(cell_type = gsub("_CIBERSORT","",cell_type),
                cell_type = gsub("_"," ",cell_type))

p1 <- cibersort_long %>% 
  ggplot(aes(ID,fraction))+
  geom_bar(stat = "identity",position = "stack",aes(fill=cell_type))+
  labs(x=NULL)+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = palette4,name=NULL)+ # iobr还给大家准备了几个色盘，贴心！
  theme_bw()+
  coord_flip()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom"
  )
p1

cibersort_long$group <- ifelse(grepl("pre$", cibersort_long$ID),
                           "pre", "post")

p2 <- ggplot(cibersort_long,aes(fct_reorder(cell_type,fraction),fraction,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  scale_fill_manual(values = palette1[c(2,4)])+ 
  theme_bw() + 
  labs(x = NULL, y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=45,hjust = 1),
        axis.text = element_text(color = "black",size = 12))+
  stat_compare_means(aes(group = group,label = ..p.format..),
                     method = "wilcox.test",label.y = 0.4)
p2


#xcell#
xcell <- deconvo_tme(eset= mRNA,
                         method = "xcell",
                         arrays = F
                         )
colnames(xcell) <- gsub("_xCell","",colnames(xcell))
colnames(xcell) <- gsub("_"," ",colnames(xcell))

xxcell <- xcell
xcell <- xcell[-c(1,2,7,13,14,20,14),]

xcell <- xcell[,c(1,5,7:15,21,22,26,28,33:35,39,50,51,55,63:65)]
xcell_long <- melt(xcell,id.vars = "ID")


cb_palette <- c("#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e",
                "#4aef7b", "#e86502", "#9ed84e", "#39ba30", "#6ad157", "#8249aa", "#99db27", "#e07233", "#ff523f",
                "#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
                "#1a918f", "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
                "#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
                "#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5","#925bea", "#63ff4f")
#也可以用iobr自带的palette 1-4

#堆叠柱状图---怎么倒过来？

ggplot(xcell_long, aes(x = ID, y = value, fill = variable)) +
  geom_bar(stat = "identity",position = "fill") +
  theme_classic() +
  coord_flip() +
  scale_fill_manual(values = palette1 ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#箱型图---对比
xcell_long$group <- ifelse(grepl("pre$", xcell_long$ID),
                           "pre", "post")

ggplot(xcell_long,aes(fct_reorder(variable,value),value,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  scale_fill_manual(values = palette1[c(2,4)])+ 
  theme_bw() + 
  labs(x = NULL, y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=90,hjust = 1),
        axis.text = element_text(color = "black",size = 12))+
  stat_compare_means(aes(group = group,label = ..p.signif..),
                     method = "kruskal.test",label.y = 0.4)

#PCA主成分----------------------------------------------------------------------
#需要mRNA，样本&基因表达的形式
pca <- read.csv("ALL_count.csv", header = T, row.names = 1)
pca <- pca[,c(1,3,8,10,12,13)]
pca <- pca[,c(1,2,3,4,6,5)]
pca <- as.data.frame(t(pca))
fenzu <- read.csv("fenzu.csv",row.names = 1)
res.pca <- prcomp(pca)
res.pca
fviz_screeplot(res.pca, addlables = TRUE)
fviz_pca_ind(res.pca, lable="none", habillage = fenzu$grp,
             addEllipses = TRUE, ellipse.level=0.95)
fviz_pca_ind(res.pca, geom.ind = "point",
             col.ind = fenzu$grp,
             addEllipses = TRUE, ellipse.type = "convex",#confidence--圆形，convex--矩形
             legend.title="Groups")


ggg <- pca[c(12,13),]
ggggg <- read.csv("ggggg.csv", row.names = 1)
rg <- prcomp(ggg)
rg
fviz_pca_ind(rg, lable="none", habillage = ggggg$grp,
             addEllipses = TRUE, ellipse.level=0.95)

#基因热图,需要先得到mRNA--------------------------------------------------------
#PD1&aAPCs，要先有mRNA这个表达矩阵
mRNA <- mRNA_count
mRNA <- mRNA[,1:6]#只要样本，基因表达量
scoregene <- readxl::read_xlsx("Score.xlsx")

pathgene <- scoregene$`Activation:Effector function`
pathgene <- na.omit(pathgene)
pathgene_gene <- mRNA[pathgene,]
pathgene_gene <- pathgene_gene[which(rowSums(pathgene_gene)>0),]#校准时一行的值不能全相同

#标准化，对照组归为1的热图：
write.csv(pathgene_gene,"3.ACTIVATION.csv")
heatmap_pathgene <- readxl::read_xlsx("4.ACTIVATION.xlsx",sheet = 2)
heatmap_pathgene <- as.data.frame(heatmap_pathgene)
rownames(heatmap_pathgene) <- heatmap_pathgene$...1
heatmap_pathgene <- heatmap_pathgene[,-1]
pheatmap(heatmap_pathgene,
         scale = "row",
         color = colorRampPalette(colors = c("#177cb0","white","#ff3300"))(100),
         fontsize = 10,
         cluster_cols = F,
         cluster_rows = F,
         angle_col = 90,border=F,
         cutree_cols = 2,
         cellwidth = 18,
         cellheight = 12,
         cutree_rows = 1)

#正常热图
pheatmap(pathgene_gene,
         scale = "row",
         fontsize = 15,
         cluster_cols = F,
         cluster_rows = T,
         angle_col = 45,border=F,
         cutree_cols = 2,
         cutree_rows = 1)

interestgene <- c("PDCD1","HAVCR2","IL17A","TIGIT","LAG3","MR1",
                  "CCR6","GZMB","GZMA","GZMH",
                  "IFNG","TNF")
interestRNA <- mRNA[interestgene,]
pheatmap(interestRNA,
         scale = "row",
         cluster_cols = F,
         cluster_rows = F,
         angle_col = 45,border=F,
         gaps_row = c(3,7))

cellcycle <- read.csv("cell cycle rpkm.csv")#还有CKR,SASP,PHENOTYPE等
cellcycle_gene <- cellcycle$SYMBOL
cellcycle_gene <- mRNA[cellcycle_gene,]
cellcycle_gene <- cellcycle_gene[which(rowSums(cellcycle_gene)>0),]
pheatmap(cellcycle_gene,
         scale = "row",
         cluster_cols = F,
         cluster_rows = T,
         angle_col = 45,border=F,
         cutree_cols = 2,
         cutree_rows = 1)

#带有标记的热图
time <- ifelse(grepl("pre$", colnames(mRNA)),
                             "pre", "post")
PPR <- ifelse(grepl("^P01|P10|P22|P18|P25", colnames(mRNA)),
              "PPR+", "nonPPR+")

MPR <- ifelse(grepl("^P01|P22|P18", colnames(mRNA)),
              "MPR", "nonMPR")
col_anno <- data.frame(row.names = colnames(mRNA),
                   PPRpositive=PPR,MPR=MPR,treatment=time)
###
PR <- ifelse(grepl("CL200106789|CL100113194|CL100113272", colnames(mRNA)),
             "PR", "SD")

col_anno <- data.frame(row.names = colnames(mRNA),
                       RECIST=PR)



interestgene <- c("CD163","CD27","GNLY","FOXP3","LAG3","CXCR6",
                  "CTLA4","CD274","PDCD1","GZMH",
                  "GZMA","GZMB","CCL5","NKG7","CXCL13","TIGIT",
                  "CD8A","PRF1","IFI16","CD4","HAVCR2","GZMK","CX3CR1","IFI44","IFI127")
interestRNA <- mRNA[interestgene,]
interestRNA <- na.omit(interestRNA)
interestRNA <- log(interestRNA,5)
interestRNA <- interestRNA+1
interestRNA[interestRNA==-Inf] <-0 
pheatmap(interestRNA,
         scale = "row",
         color = colorRampPalette(colors = c("#177cb0","white","#ff3300"))(100),
         annotation_col =col_anno,
         fontsize = 10,
         cluster_cols = T,
         cluster_rows = T,
         angle_col = 90,border=F,
         cutree_cols = 2,
         cellwidth = 12,
         cellheight = 12,
         cutree_rows = 1)

#基因表达箱型图-----------------------------------------------------------------
gene <- c("CD247","CTLA4","GZMB","NKG7","PRF1","OASL","IFI44")

gene <- c("CD8A","GZMH","GZMA","NKG7","CD274","FOXP3","IDO1","CD163")

gene <- c("CXCL9","CXCL10","CXCL11","FAS")

gene <- c("MMP1","OLFM4","IFI27","PDZK1IP1","LCN2","S100A9")
gene <- as.vector(gene)
mRNA <- log(mRNA,5)
mRNA <- mRNA+1
boxgene <- mRNA[gene,]
boxgene <- t(boxgene)
boxgene <- as.data.frame(boxgene)
Exp_plot <- boxgene

###
PR <- ifelse(grepl("CL200106789|CL100113194|CL100113272", colnames(mRNA)),
             "PR", "SD")

col_anno <- data.frame(row.names = colnames(mRNA),
                       RECIST=PR)

col_anno$sample <- rownames(col_anno)
Exp_plot<- Exp_plot[col_anno$sample,]

Exp_plot$sam=col_anno$RECIST
Exp_plot$sam <- factor(Exp_plot$sam,levels=c("PR","SD"))

#循环语句实现一图多出
col <-c("#5CB85C","#337AB7","#F0AD4E","#D9534F")
plist2<-list()
for (i in 1:length(gene)){
  bar_tmp<-Exp_plot[,c(gene[i],"sam")]
  colnames(bar_tmp)<-c("Expression","sam")
  my_comparisons1 <- list(c("PR", "SD")) 
  pb1<-ggboxplot(bar_tmp,
                 x="sam",
                 y="Expression",
                 color="sam",
                 bxp.errorbar.width = 2,
                 width = 0.5,
                 size=1,
                 font.label = list(size=30), 
                 palette = col)+theme(panel.background =element_blank())
  pb1<-pb1+theme(axis.line=element_line(colour="black"))+theme(axis.title.x = element_blank())
  pb1<-pb1+theme(axis.title.y = element_blank())+theme(axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1))
  pb1<-pb1+theme(axis.text.y = element_text(size = 15))+ggtitle(gene[i])+theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))
  pb1<-pb1+theme(legend.position = "NA")#
  pb1<-pb1+stat_compare_means(hide.ns = F,
                              comparisons =c(my_comparisons1),
                              label="p.signif")
  plist2[[i]]<-pb1 
}
plot_grid(plist2[[1]],plist2[[2]],plist2[[3]],
          plist2[[4]],plist2[[5]],plist2[[6]],ncol = 3)
