
library(BiocManager)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(gridExtra)
library(readxl)
#library(Matrix)
library(pheatmap)
library(ggrepel)
library(ggsignif)
library(ggpubr)
library(ggridges)
library(harmony)
library(monocle)
library(Biobase)


data <- read.table("GSE120575_patient_ID_single_cells.txt",header = T,sep = "\t",skip = 19)
Exp <- as.data.frame(data.table::fread(file  = "GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt"))

#Standared 10x genomic dataset
pbmc.data <- Read10X(data.dir = "~/Rdata/single_cell_SEQ/filtered_gene_bc_matrices/hg19/")
example <- CreateSeuratObject(counts = pbmc.data, project = "pmbc3k", min.cells = 3, min.features = 200)
#过滤检测少于200个基因的细胞（min.features = 200）和少于3个细胞检测出的基因（min.cells = 3）

#read--matrix.gz
NSCLC <- read.table(gzfile("GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"), header = FALSE, row.names = 1, sep = "\t")
metadata <- read_excel("GSE207422_NSCLC_scRNAseq_metadata.xlsx")
metadata <- metadata[1:15,]
colnames(NSCLC) <- NSCLC[1,]

select <- NSCLC[,grepl("^BD_immune(01|02|03|04|05|08)",colnames(NSCLC))]#grepl!!
rm(NSCLC)

write.csv(select,file = "select.csv")
select <- read.csv(file = "select.csv")
colnames(select) <- select[1,]
select <- select[-1,]
rownames(select) <- select[,1]
select <- select[,-1]
select[] <- as.data.frame(lapply(select,as.numeric))

pbmc <- CreateSeuratObject(counts = select,project = "GSE207422",
                           min.cells = 3, min.features = 200)


#add metadata according to the meta-infomation
pbmc <- AddMetaData(object = pbmc,
                             metadata = ifelse(grepl("^BD_immune01|BD_immune05|BD_immune08", colnames(pbmc)),
                                               "pre", "post"), col.name = "time")


TEST <- AddMetaData(object = TEST, metadata = ifelse(grepl("^P", TEST$orig.ident), TEST$orig.ident, 
                      ifelse(TEST$orig.ident == "CL100100202", "P14pre",
                      ifelse(TEST$orig.ident == "CL100113194", "P12post",
                      ifelse(TEST$orig.ident == "CL100113271", "P16pre",
                      ifelse(TEST$orig.ident == "CL100113272", "P13post",
                      ifelse(TEST$orig.ident == "CL100120780", "P15pre",
                      ifelse(TEST$orig.ident == "CL100122110", "P15post",
                      ifelse(TEST$orig.ident == "CL200106789", "P11post",
                      ifelse(TEST$orig.ident == "CL200124680", "P16post", NA))))))))), col.name = "patient")



TEST <- AddMetaData(object = TEST,
                    metadata = ifelse(grepl("pre", TEST@meta.data$patient),
                                      "pre", "post"), col.name = "time")
#先取治疗后的再比较
TEST <- AddMetaData(object = TEST,
                    metadata = ifelse(grepl("P11post|P12post|P13post", TEST@meta.data$patient),
                                      "PR", "SD"), col.name = "RECIST")


# 获取细胞名字中的数字
cell_num <- as.numeric(gsub("BD_immune(\\d+)_.*", "\\1", colnames(TEST)))
# 将数字转换为patient编号
patient <- sprintf("patient%02d", cell_num)
# 添加metadata
TEST <- AddMetaData(TEST, metadata = patient, col.name = "patient")


saveRDS(pbmc,file = "SELECT.rds")
save.image(file = "data_harmony.RData")

#清空环境
rm(list = ls())
gc()

#质量控制-----------------------------------------------------------------------

        ###低质量/濒死细胞常表现出广泛的线粒体污染
        ###使用PercentageFeatureSet()函数计算线粒体QC指标

#对TEST群重新降维,循环----------------------------------------------------------
TEST[["percent.mt"]] <- PercentageFeatureSet(TEST, pattern = "^MT-")
#展示前5个细胞的QC指标
head(TEST@meta.data, 5)
#使用小提琴图可视化QC指标(基因表达数/.基因表达量之和/线粒体基因占比)
VlnPlot(TEST, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#看nC和nF以及mt之间的相关性
plot1 <- FeatureScatter(TEST, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(TEST, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#过滤线粒体基因比值过高的细胞
TEST <- subset(TEST, subset = nCount_RNA > 200 & nCount_RNA < 30000 & percent.mt < 180)#!!!!!!!!!!!!!!

#数据标准化
TEST <- NormalizeData(TEST, normalization.method = "LogNormalize", scale.factor = 10000)#!!!!!!!!!!!!!!

#鉴定高细胞间变异的特征基因
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 5000)
# 查看最高变的10个基因
top10 <- head(VariableFeatures(TEST), 10)
plot1 <- VariableFeaturePlot(TEST,pt.size = 0.5, selection.method = "vst")
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
plot2

all.genes <- rownames(TEST)
TEST <- ScaleData(TEST, features = all.genes)

#线性降维，PCA主成分分析
TEST <- RunPCA(TEST, features = VariableFeatures(object = TEST),verbose = T)
TEST  <-  RunHarmony(TEST,group.by.vars = "patient", plot_convergence = TRUE, lambda= 1.4)#harmony 去除批次效应
DimHeatmap(TEST, dims = 1:20, cells = 500, balanced = TRUE)

ElbowPlot(TEST,ndims = 30)#看拐点，ndim= 20-50
#print(TEST[["pca"]], dims = 1:5, nfeatures = 5)
#VizDimLoadings(TEST, dims = 1:2, reduction = "pca")
#JackStrawPlot(TEST, dims = 1:15)
#非线性降维（UMAP/tSNE）
#UMAP
TEST <- RunUMAP(TEST, dims = 1:10)#,reduction = "harmony"
#细胞聚类
TEST <- FindNeighbors(TEST,dims = 1:10)#,reduction = "harmony"
#使用真实信号的前6个PC分类  #TEST <- FindNeighbors(TEST, k.param = 30)
TEST <- FindClusters(TEST, resolution = 0.6)#resolution3k个细胞时0.4-1.2，越大分类越多
DimPlot(TEST,reduction = 'umap',group.by = 'seurat_clusters',label = F,raster=FALSE)
DimPlot(TEST,reduction = 'umap',group.by = 'RECIST',label = F)#split.by:post,pre!!!!!!!!!!!
DimPlot(TEST,reduction = 'umap',split.by = 'orig.ident',label = F)
DimPlot(TEST,reduction = 'tsne',group.by = 'celltype',label = T,raster=FALSE)
dif <- TEST[,TEST@meta.data$seurat_clusters %in%
              c(5,17)]
FeaturePlot(TEST,features = "cnvscore",reduction = 'umap',
            cols = rev(brewer.pal(6, name = "GnBu")),order = T,pt.size = 1)

marker <- c("CLDN4","CLDN7")
FeaturePlot(TEST,features = "",reduction = 'umap')

VlnPlot(TEST,features = marker,group.by = "seurat_clusters")

marker <- c("LGALS1","LGALS3","LGALS9")
DotPlot(TEST,features = marker)

TEST <- TEST[,TEST@meta.data$seurat_clusters %in%
              c(0,1,2,4,5)]
x8 <- TEST[,TEST@meta.data$seurat_clusters %in%
             c(1)]
saveRDS(x8,file = "final.rds")

#循环#########################################----------------------------------
rm(CD8)
TEST <- MAIT
rm(MAIT)
FeaturePlot(post, features = c("GZMK", "CD8A", "CD8B"))#CD8+T
marker <- c("GZMK", "CD8A", "CD8B")
DotPlot(TEST,features = marker)

saveRDS(TEST, file = "filtered.rds")
save.image(file = "data.RData")
#展示细胞比例-------------------------------------------------------------------
#首先取子集分开两者
pre <- TEST[,TEST@meta.data$time %in%
              c("pre")]
pre <- RunUMAP(pre, dims = 1:6)
DimPlot(pre,reduction = 'umap',group.by = 'seurat_clusters',label = T)
maitpre <- WhichCells(pre,expression = SLC4A10>0)
ll <- pre[,colnames(pre) %in%
             maitpre]
#maitpre <- subset(pre,pre[["RNA"]]@data["SLC4A10",]>0)


#已经注释了细胞类型，看样本中比例情况
Cellratio <- prop.table(table(Idents(TEST),TEST$patient),margin = 2)
Cellratio <- as.data.frame(Cellratio)

library(reshape2)
library(IOBR)
library(tidyverse)
cellper <- dcast(Cellratio,Var2~Var1, value.var = "Freq")#长数据转为宽数据
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]
cellper$ID <- rownames(cellper)
cell_long <- melt(cellper,id.vars = "ID")

cell_long$group <- ifelse(grepl("post", cell_long$ID),
                          "post", "pre")
cell_long$group <- ifelse(grepl("P11post|P12post|P13post", cell_long$ID),
                          "PR", "SD")

ggplot(cell_long,aes(fct_reorder(variable,value),value,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  scale_fill_manual(values = palette1[c(2,4)])+ 
  theme_bw() + 
  labs(x = NULL, y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=90,hjust = 1),
        axis.text = element_text(color = "black",size = 12))+
  stat_compare_means(aes(group = group,label = ..p.format..),
                     method = "kruskal.test",label.y = 0.4)



#注释细胞类型（根据经验）-------------------------------------------------------
maker <- c("CD3D","CD2","IL7R","BATF","IL2RA","FOXP3","GZMA","CD8A",
           "CD8B","UBE2C","TOP2A","CDK1","MS4A1","CD79A","CD19","MZB1","DERL3",
           "SDC1","CD68","C1QC","S100A9","S100A8","FCN1","CD1C","CLEC10A","CCL22",
           "LAMP3","TPSB2","VWF","PECAM1","DCN","LUM","COL6A3","BGN","MYL9","ACTA2")

maker <- c("EPCAM","KRT19","KRT18","CD68","CSF1R","CD14","CD2","CD3D","CD3E",
           "CD79A","CD19","MS4A1","CLDN5","PECAM1","VWF","TPSAB1","TPSB2","KIT",
           "COL1A1","COL3A1","DCN")
DotPlot(TEST, features = maker ,
        assay='RNA' ) +  #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))

new.cluster.ids <- c("T cell", "T cell", "B cell", "epithelial cell", "epithelial cell", "mast cell",
                     "fibroblast","myeloid cell","endothelial cell","fibroblast","epithelial cell",
                     "epithelial cell","endothelial cell")



new.cluster.ids <- c("malignant", "normal", "malignant-like", "malignant", "malignant",
                     "normal-like", "malignant-like", "malignant", "malignant-like", "malignant",
                     "malignant", "normal-like", "malignant-like", "malignant", "malignant",
                     "malignant-like", "malignant", "malignant-like", "malignant", "malignant")  

names(new.cluster.ids) <- levels(TEST)
TEST <- RenameIdents(TEST, new.cluster.ids)

TEST$celltype <- Idents(TEST)

#基因特征（耗竭/活化等）打分函数！！--------------------------------------------
scoregene <- readxl::read_xlsx("~/wyc_20231021/single_cell_SEQ/project/NSCLC/Score.xlsx")
library(ggsci)
library(cowplot)
path <- c("Activation:Effector function","Chemokine/Chemokine receptor","Cytotoxicity","Exhaustion")
plist <- list()
for (i in 1:length(path)){
pathgene <- scoregene[[path[i]]]
pathgene <- as.data.frame(pathgene)
pathgene <- na.omit(pathgene)
pathgene <- as.list(pathgene)
TEST <- Seurat::AddModuleScore(TEST,features = pathgene,name = paste0(gsub("[^[:alnum:]]", ".", path[i])))
boxdata <- FetchData(TEST,vars = c("time",paste0(gsub("[^[:alnum:]]", ".", path[i]), "1")))
my_comparisons <- list( c("pre", "post"))
P <- ggplot(boxdata,aes_string("time",paste0(gsub("[^[:alnum:]]", ".", path[i]), "1")))
P <- P+geom_boxplot(aes(fill = time))+theme_bw()+RotatedAxis()+
  geom_signif(comparisons = my_comparisons,
              step_increase = 0.3,
              map_signif_level = F,
              test = wilcox.test)+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  scale_fill_npg()+
  theme(axis.text.y = element_text(size = rel(1),colour = "black"),
        axis.text.x = element_text(size = rel(1),colour = "black",angle = 0),
        legend.title = element_text(size = rel(1),colour = "black"),
        legend.text = element_text(size = rel(1),colour = "black"),
        axis.title = element_text(size = rel(1),colour = "black"))+
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test")
plist[[i]] <- P
}
plot_grid(plist[[1]],plist[[2]],plist[[3]],
          plist[[4]],ncol = 2)#,plist[[5]],plist[[6]]

#MAIT
MAITgene <- list(c("TRAV1-2","SLC4A10","KLRB1","RORC","ZBTB16"))

taggene <- list(c("CLDN4","CLDN7","TFF3","REG4"))
TEST <- Seurat::AddModuleScore(TEST,features = taggene,name = "tag_score")

FeaturePlot(TEST,  features ="tag_score1" ,reduction = "umap",
            cols = rev(brewer.pal(10, name = "RdBu")),order = T,pt.size = 1)

VlnPlot(TEST,group.by = "time",features = "MAIT_score1")

boxdata <- FetchData(TEST,vars = c("RECIST","MAIT_score1"))
boxdata <- boxdata[boxdata$MAIT_score1>0,]
P <- ggplot(boxdata,aes(RECIST,MAIT_score1))
my_comparisons <- list( c("PR", "SD"))
P+geom_boxplot()+theme_bw()+RotatedAxis()+
  geom_signif(comparisons = my_comparisons,
              step_increase = 0.3,
              map_signif_level = F,
              test = t.test)+
  theme(axis.text.y = element_text(size = rel(1.5),colour = "black"),
        axis.text.x = element_text(size = rel(3),colour = "black"),
        legend.title = element_text(size = rel(1.5),colour = "black"),
        legend.text = element_text(size = rel(1.5),colour = "black"),
        axis.title = element_text(size = rel(1.5),colour = "black"))+
  stat_compare_means(comparisons = my_comparisons,
                  method = "kruskal.test")


#heatmap只需要基因向量即可------------------------------------------------------
DoHeatmap(epicell,group.by = 'celltype',maker,#pathgene must be chr frame
          assay = "RNA")+scale_fill_gradientn(colors = c("navy","white","firebrick3"))

#get counts first!然后可以用heatmap实现聚类-------------------------------------
mRNA <- GetAssayData(object = TEST, slot = "scale.data")
mRNA <- as.data.frame(mRNA)

top_30 <- mRNA[pathgene,]
top_30 <- top_30[which(rowSums(top_30)>0),]

pheatmap(top_30,
         scale = "row",
         color = colorRampPalette(colors = c("#177cb0","white","#ff3300"))(100),
         fontsize = 10,
         cluster_cols = T,
         cluster_rows = T,
         angle_col = 90,border=F,
         cutree_cols = 2,
         cellwidth = 0.3,
         cellheight = 10,
         cutree_rows = 1)

#细胞亚群分类后差异表达分析-----------------------------------------------------

#某一群与其他的比较
DEG <- FindMarkers(TEST,ident.1 = "post", ident.2 = "pre",group.by = "time") 
write.csv(DEG,file = "DEG_MAIT_time.csv")
DEG_RS <- FindMarkers(TEST,ident.1 = "PR", ident.2 = "SD",group.by = "RECIST")
write.csv(DEG_RS,file = "DEG_RS_result.csv")

my_result <- read.csv(file = "DEG_RS_result.csv")

colnames(my_result) <- c("Gene_symbol","pvalue","log2FoldChange","pct.1","pct.2","padj")

my_result$regulate <- ifelse(my_result$pvalue > 0.05, "unchanged",
                             ifelse(my_result$log2FoldChange > 0.5, "up-regulated",
                                    ifelse(my_result$log2FoldChange < -0.5, "down-regulated","unchanged")))
table(my_result$regulate)

DEG_deseq2 <- subset(my_result, pvalue < 0.1 & abs(log2FoldChange) > 1)#######注意要改
upgene <- DEG_deseq2[DEG_deseq2$regulate=='up-regulated',]
downgene <- DEG_deseq2[DEG_deseq2$regulate=='down-regulated',]

# 找出每个cluster的标记与所有剩余的细胞相比较，只报告阳性细胞
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% slice_max(n=2,order_by = avg_log2FC)
#每个聚类前10个差异基因表达热图(如果小于10，则绘制所有标记)
top10_cluster <- MAIT %>% top_n(n = 30, wt = avg_log2FC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

deg_top30 <- rownames(DEG[1:30,])

pathgene <- deg_top30

deg_top30 <- as.data.frame(deg_top30)
rownames(deg_top30) <- deg_top30$deg_top30

#火山图-------------------------------------------------------------------------
plot(my_result$log2FoldChange,-log2(my_result$pvalue))

###可以改一下regulate的cutoff条件
my_result$regulate <- ifelse(my_result$pvalue > 0.2, "unchanged",
                             ifelse(my_result$log2FoldChange > 0.2, "up-regulated",
                                    ifelse(my_result$log2FoldChange < -0.2, "down-regulated","unchanged")))

PWYL <- scoregene$Exhaustion
PWYL <- na.omit(PWYL)
PWYL_gene <- my_result[PWYL,]

PWYL_gene$Gene_symbol <- row.names(PWYL_gene)
PWYL_gene <- na.omit(PWYL_gene)

#title！show指定gene/show界限gene
P <- ggplot(data=PWYL_gene, aes(x = log2FoldChange, y = -log10(pvalue), color=regulate)) +
  geom_point(shape=16, size=3) +
  theme_set(theme_set(theme_bw(base_size = 20))) +
  xlab("log2 FoldChange") +
  ylab("-log10 p-value") +
  theme(plot.title = element_text(size = 10, hjust = 2.5)) +
  theme_classic() +
  scale_color_manual(values = c('navy','gray','firebrick')) +  #'#177cb0','gray','#ff3300'
  geom_vline(xintercept = c(-0.4,0.4), lty=4, col='gray', lwd=0.8) +#logFC分界???
  geom_hline(yintercept = -log10(0.1), lty=2, col='gray', lwd=0.6) +#adj.p.val分界???
  labs(title = 'cytotoxicity') +
  geom_text_repel(data = filter(PWYL_gene),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 10),
                  aes(label = Gene_symbol, ),min.segment.length = 0,force = 100,
                  size = 5)+
  coord_cartesian(xlim=c(-max(abs(PWYL_gene$log2FoldChange)), max(abs(PWYL_gene$log2FoldChange))))
plot(P)

#show all DEgene
P <- ggplot(data=my_result, aes(x = log2FoldChange, y = -log10(pvalue), color=regulate)) +
  geom_point(shape=16, size=3) +
  theme_set(theme_set(theme_bw(base_size = 20))) +
  xlab("log2 FoldChange") +
  ylab("-log10 p-value") +
  theme(plot.title = element_text(size = 10, hjust = 2.5)) +
  theme_classic() +
  scale_color_manual(values = c('#86D47D','#DEF0DF','#C43B99')) +
  geom_vline(xintercept = c(-1,1), lty=4, col='gray', lwd=0.8) +#logFC分界???
  geom_hline(yintercept = 20, lty=2, col='gray', lwd=0.6) +#adj.p.val分界???-log10()
  labs(title = 'PD1-aAPCs') +
  geom_text_repel(data = filter(my_result, abs(my_result$log2FoldChange) > 1.5 & -log10(my_result$pvalue) > 20),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 10),
                  aes(label = Gene_symbol, ),min.segment.length = 0,force = 1,
                  size = 5) 
plot(P)


#取子集-------------------------------------------------------------------------
pre <- pbmc[,pbmc@meta.data$time %in%
                      c("Pre")]
TEST <- TEST[,TEST@meta.data$patient %in%
               c("P9post","P11post","P12post","P13post","P14post","P15post","P16post")]

TEST <- subset(TEST,subset = SLC4A10>0 )

pre <- RunUMAP(pre, dims = 1:6)
FeaturePlot(pre,features = c("SLC4A10"),reduction = 'umap')
DimPlot(pre, reduction = "umap", label = TRUE)

TEST <- TEST[,TEST@meta.data$celltype %in%
               c("CD8+T cell")]

epicell <- TEST[,TEST@meta.data$celltype %in%
               c("Epithelial cell")]


tcell <- TEST[,TEST@meta.data$celltype %in%
                c("CD8+T cell","CD4+T cell")]
bcell <- TEST[,TEST@meta.data$celltype %in%
                c("B cell")]
epithelialcell <- TEST[,TEST@meta.data$celltype %in%
                c("epithelial cell")]
myeloidcell <- TEST[,TEST@meta.data$celltype %in%
                c("myeloid cell")]
mastcell <- TEST[,TEST@meta.data$celltype %in%
                c("mast cell")]
endothelialcell <- TEST[,TEST@meta.data$celltype %in%
                c("endothelial cell")]
fibroblast <- TEST[,TEST@meta.data$celltype %in%
                c("fibroblast")]

TEST <- tcell
TEST <- bcell
TEST <- epithelialcell
TEST <- myeloidcell
TEST <- mastcell
TEST <- endothelialcell
TEST <- fibroblast


TEST <- merge(bcell,y = c(endothelialcell,epithelialcell,fibroblast,
                          mastcell,myeloidcell,tcell),
              add.cell.ids = c("bcell","endothelialcell","epithelialcell",
                               "fibroblast","mastcell","myeloidcell","tcell"),
              project = "combined")




#cellchat-----------------------------------------------------------------------
library(CellChat)
library(Seurat)
library(dplyr)
library(SeuratData)
library(patchwork)
library(ggplot2)
library(ggalluvial)
library(svglite)
options(stringsAsFactors = F)

CellChatDB <- CellChatDB.human
CellChatDB$interaction[1:4,1:4]
unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")#此外还有ECM（细胞外基质），cellcell-contact


###第一种

PR <- TEST[,TEST@meta.data$RECIST %in%
               c("PR")]
PR_chat <- createCellChat(object=PR,meta= PR@meta.data,group.by = "celltype")
PR_chat
summary(PR_chat)
str(PR_chat)
levels(PR_chat@idents)
PR_groupSize <- as.numeric(table(PR_chat@idents))
PR_chat@DB <- CellChatDB.use

## 在矩阵的所有的基因中提取signaling gene，结果保存在data.signaling(13714个基因，过滤完只有270个）
PR_chat <- subsetData(PR_chat)
future::plan("multiprocess", workers = 4)

PR_chat <- identifyOverExpressedGenes(PR_chat)#相当于Seurat的FindMarkers，找每个细胞群中高表达的配体受体
PR_chat <- identifyOverExpressedInteractions(PR_chat) 
#上一步运行的结果储存在cellchat@LR$LRsig
PR_chat <- projectData(PR_chat, PPI.human) 
#找到配体受体关系后，projectData将配体受体对的表达值投射到PPI上，来对@data.signaling中的表达值进行校正。结果保存在@data.project
#根据表达值推测细胞互作的概率（cellphonedb是用平均表达值代表互作强度）。
PR_chat <- computeCommunProb(PR_chat) #如果不想用上一步PPI矫正的结果，raw.use = TRUE即可。, raw.use = TRUE, population.size = TRUE
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
PR_chat <- filterCommunication(PR_chat, min.cells = 10)
df.net <- subsetCommunication(PR_chat)
write.csv(df.net, "net_PR.csv")
PR_chat <- computeCommunProbPathway(PR_chat)
df.netp <- subsetCommunication(PR_chat, slot.name = "netP")
write.csv(df.netp, "net_PR_pathway.csv")
PR_chat <- aggregateNet(PR_chat)

groupSize <- as.numeric(table(PR_chat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(PR_chat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions")
netVisual_circle(PR_chat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength")

###第二种
SD <- TEST[,TEST@meta.data$RECIST %in%
             c("SD")]
SD_chat <- createCellChat(object=SD,meta= SD@meta.data,group.by = "celltype")
SD_chat
summary(SD_chat)
str(SD_chat)
levels(SD_chat@idents)
SD_groupSize <- as.numeric(table(SD_chat@idents))
SD_chat@DB <- CellChatDB.use

## 在矩阵的所有的基因中提取signaling gene，结果保存在data.signaling(13714个基因，过滤完只有270个）
SD_chat <- subsetData(SD_chat)
future::plan("multiprocess", workers = 4)

SD_chat <- identifyOverExpressedGenes(SD_chat)#相当于Seurat的FindMarkers，找每个细胞群中高表达的配体受体
SD_chat <- identifyOverExpressedInteractions(SD_chat) 
#上一步运行的结果储存在cellchat@LR$LRsig
SD_chat <- projectData(SD_chat, PPI.human) 
#找到配体受体关系后，projectData将配体受体对的表达值投射到PPI上，来对@data.signaling中的表达值进行校正。结果保存在@data.project
#根据表达值推测细胞互作的概率（cellphonedb是用平均表达值代表互作强度）。
SD_chat <- computeCommunProb(SD_chat) #如果不想用上一步PPI矫正的结果，raw.use = TRUE即可。, raw.use = TRUE, population.size = TRUE
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
SD_chat <- filterCommunication(SD_chat, min.cells = 10)
df.net <- subsetCommunication(SD_chat)
write.csv(df.net, "net_SD.csv")
SD_chat <- computeCommunProbPathway(SD_chat)
df.netp <- subsetCommunication(SD_chat, slot.name = "netP")
write.csv(df.netp, "net_SD_pathway.csv")
SD_chat <- aggregateNet(SD_chat)

groupSize <- as.numeric(table(SD_chat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(SD_chat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions")
netVisual_circle(SD_chat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength")



object.list <- list(PR = PR_chat, SD = SD_chat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

#对比二者之间的不同的交互数量及权重
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

#将二者对比，以circle形式展示细胞种类间交互的变化，前者比后者（红色增强，蓝色减弱）
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

#热图形式展示，左侧是发出者，下面是接收者###
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

#二维气泡图###
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]], slot.name = "netP")
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)

netVisual_bubble(cellchat, sources.use =c(1:9), targets.use = c(5),  comparison = c(1, 2), angle.x = 90)


#在某一个信号通路中比较细胞群间的区别
object.list[[1]]@netP$pathways  #查看都有哪些信号通路

pathways.show <- c("CXCL") 

#弦图
par(mfrow = c(1,1), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}


par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))


#比较每个细胞群传入传出信号
library(ComplexHeatmap)

i = 1
# combining all the identified signaling pathways from different datasets
#OUTGOING
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
#INCOMING
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
#OVERALL
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

saveRDS(cellchat,file = "all_cellchat.rds")
saveRDS(object.list,file = "all_cellchat_objectlist.rds")



