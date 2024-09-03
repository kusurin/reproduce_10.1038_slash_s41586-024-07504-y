seu_obj <- Seurat::CreateSeuratObject(Seurat::Read10X("../",gene.column = 1)) #or subset
seu_obj <- Seurat::CreateSeuratObject(Seurat::Read10X("../",gene.column = 1),meta.data = data.frame(library = seu_meta$library,row.names = seu_meta$barcode))

#筛低质量组：feature和count和线粒体基因与其他比起来有显著差异，snRNA筛掉了2个

seu_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu_obj,pattern = "^mt-") #全是0？淦
print(grep("^MT-",rownames(seu_obj),value=TRUE))
print(grep("^MT",rownames(seu_obj),value=TRUE))
print(grep("^mt",rownames(seu_obj),value=TRUE))

#应该先ANOVA吧
anova_oridata <- seu_obj@meta.data[, c("nCount_RNA", "nFeature_RNA", "percent.mt", "orig.ident")]

for(ident in unique(anova_oridata$orig.ident)) {
  qqnorm(anova_oridata[anova_oridata$orig.ident == ident,]$percent.mt,main=paste("(nor)QQ-mt-",ident))
  qqline(anova_oridata[anova_oridata$orig.ident == ident,]$percent.mt)
}

install.packages("car")

car::leveneTest(percent.mt ~ orig.ident,data=anova_oridata)

#绷

anova_count <- aov(nCount_RNA ~ orig.ident,data = anova_oridata)
summary(anova_count)

kruskal.test(nCount_RNA ~ orig.ident,data=anova_oridata)

oneway.test(nCount_RNA ~ orig.ident,data=anova_oridata,var.equal = FALSE)

install.packages("PMCMRplus")

GH_count <- gamesHowellTest( data = anova_oridata,nCount_RNA ~ orig.ident)

means <- anova_oridata %>%
  group_by(orig.ident) %>%
  summarise(count = mean(nCount_RNA, na.rm = TRUE))

#中位数sex和age差得最多，显著性drug和sex差得最多

#不对，看了meta，这七个是大组，里面还有小组，好像应该筛小组
#刚好435099，不质控了以上全部略过（）

#标准化
seu_obj <- Seurat::NormalizeData(seu_obj,normalization.method="LogNormalize",scale.factor=10000)
#挑选高变异基因
seu_obj <- Seurat::FindVariableFeatures(seu_obj,selection.method="vst",nFeatures = 2000)
p_temp <- Seurat::VariableFeaturePlot(seu_obj)
p_temp <- Seurat::LabelPoints(plot = p_temp, points = head(Seurat::VariableFeatures(seu_obj), 10) , repel = TRUE, size=2.5) + ggplot2::theme(legend.position = c(0.1,0.8))
#PCA
seu_obj <- Seurat::ScaleData(seu_obj)
seu_obj <- Seurat::RunPCA(seu_obj, features = Seurat::VariableFeatures(seu_obj))
p_temp <- Seurat::DimPlot(seu_obj, reduction = "pca", group.by="library")
#整合
seu_obj <- RunHarmony(seu_obj, group.by.vars = "orig.ident", dims.use = 1:50) #default pca
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
seu_obj <- Seurat::RunUMAP(seu_obj,reduction="harmony",dims=1:50)

Seurat::DimPlot(seu_obj, reduction = "umap", group.by = "seurat_clusters")

library <- seu_obj@meta.data$library
main_library <- sub("_[^_]*$", "",library)
main_library <- as.factor(main_library)
seu_obj@meta.data$sub_orig.ident <- main_library
Seurat::DimPlot(seu_obj, reduction = "umap", group.by = "sub_orig.ident")

#手动注释细胞类型