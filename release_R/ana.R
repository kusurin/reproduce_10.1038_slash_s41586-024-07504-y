sn_meta <- read.table("../data_origin/GSE234774_rnaseq_meta.txt.gz",sep = "\t",header = T)
seu_obj <- Seurat::CreateSeuratObject(Seurat::Read10X("../",gene.column = 1))
seu_obj$library <- as.factor(sn_meta$library)
seu_obj$meta_layer1 <- as.factor(sn_meta$layer1)
seu_obj$meta_layer2 <- as.factor(sn_meta$layer2)


sce <- seu_obj@assays$RNA$counts
sce <- decontX::decontX(x=sce)
seu_obj$contamination_decontX <- sce$contamination
seu_obj<-subset(seu_obj,subset = nCount_RNA >= 200)
#筛低质量组：feature和count和线粒体基因与其他比起来有显著差异，snRNA筛掉了2个
seu_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu_obj,pattern = "^mt-") #全是0？淦

print(grep("^mt",rownames(seu_obj),value=TRUE,ignore.case = TRUE,perl = TRUE))
#最多82.58%
print(grep("^hb(?!p1|egf).*",rownames(seu_obj),value=TRUE,ignore.case = TRUE,perl = TRUE))
seu_obj[["percent.hb"]] <- Seurat::PercentageFeatureSet(seu_obj,pattern = "^Hb[^ep]") #不支持perl嗷

#seu_obj<-subset(seu_obj,subset = percent.mt<20 & percent.hb <3)

#标准化
seu_obj <- Seurat::NormalizeData(seu_obj,normalization.method="LogNormalize",scale.factor=10000)
#挑选高变异基因
seu_obj <- Seurat::FindVariableFeatures(seu_obj,selection.method="vst",nFeatures = 2000)
#p_temp <- Seurat::VariableFeaturePlot(seu_obj)
#p_temp <- Seurat::LabelPoints(plot = p_temp, points = head(Seurat::VariableFeatures(seu_obj), 10) , repel = TRUE, size=2.5) + ggplot2::theme(legend.position = c(0.1,0.8))
#PCA
seu_obj <- Seurat::ScaleData(seu_obj)
seu_obj <- Seurat::RunPCA(seu_obj, features = Seurat::VariableFeatures(seu_obj))
#p_temp <- Seurat::DimPlot(seu_obj, reduction = "pca", group.by="library")

#Seurat::ElbowPlot(seu_obj,ndims = 50)

#png("DimHeatmap60.png", width = 2000, height = 1500)
#Seurat::DimHeatmap(seu_obj,dims = 1:50, cells = 500, balanced = TRUE,ncol = 10,raster=FALSE)
#dev.off()
#整合
seu_obj <- harmony::RunHarmony(seu_obj, group.by.vars = "orig.ident", dims.use = 1:50) #default pca
#seu_obj <- Seurat::IntegrateLayers(
#  object = seu_obj, method = HarmonyIntegration,
#  orig.reduction = "pca", new.reduction = "harmony",
#  verbose = FALSE
#)

#Leiden法聚类
seu_obj <- Seurat::FindNeighbors(seu_obj,reduction="harmony",dims=1:50)
reticulate::py_install("leidenalg")
seu_obj <- Seurat::FindClusters(seu_obj, resolution = 0.05, method = "igraph",algorithm = 4) #默认直接virt1400g
#事实上，根据https://github.com/satijalab/seurat/issues/2294，Leiden法或许会将稀疏矩阵转换成稠密矩阵，可以尝试使用method = igraph

#用umap好了(确实用的也是umap)
seu_obj <- Seurat::RunUMAP(seu_obj,reduction="harmony",dims=1:50,min.dist=0.2)#这个参数好看一点


Seurat::DimPlot(seu_obj, reduction = "umap", group.by = "seurat_clusters",raster = FALSE,label = T)

library <- seu_obj@meta.data$library
main_library <- sub("_[^_]*$", "",library)
main_library <- as.factor(main_library)
seu_obj@meta.data$sub_orig.ident <- main_library
Seurat::DimPlot(seu_obj, reduction = "umap", group.by = "sub_orig.ident")
