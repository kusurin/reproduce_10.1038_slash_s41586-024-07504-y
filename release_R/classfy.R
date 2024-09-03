
#手动注释细胞类型
markers <- list(
  oligodendrocyte = c("Mobp", "Cldn11", "Plp1", "Mog", "Mag"),
  neuron = c("Meg3", "Snhg11", "Snap25", "Scg2","Stmn3"),
  astro = c("Gfap", "Aqp4", "Slc1a3", "Slc1a2"),
  ependymal = c("Foxj1", "Tmem212", "Ccdc153", "Cfap44"),
  fibroblast = c("Col3a1", "Col1a1", "Clec3b", "Dcn"),
  endothelial = c("Pecam1", "Cldn5", "Ctla2a", "Slco1a4"),
  immune = c(
    "Retnlg", "Ngp", "Stfa2l1", # Neutrophil
    "Ifi27l2a", "Crip1", "Ifitm3", # Monocyte
    "Lyz2", "Pf4", "F13a1", # Macrophage
    "Irf8", "Cd74", "Cd83", # Dendritic
    "Cd3g", "Cd3d", "Lef1", # Lymphocyte
    "Tmem119", "Tlr7", "Clec5a", # Microglia
    "Trbc2", "Trac", "Trdc", # T cell
    "Iglc1","Pax5","Fcrla"#B
  )
)

for(type in markers){      
  cat("\33[34m",type,"\33[0m\n")
  for (genename in type) {
    i<-grep(paste0("^",genename,"$"),rownames(seu_obj),value=TRUE,ignore.case = TRUE,perl = TRUE)
    if(length(i)==0){
      cat("\33[31m",genename,"not found\33[0m\n")
    }
    else print(genename)
  }
}

for(type in markers){      
  cat("\33[34m",type,"\33[0m\n")
  for (genename in type) {
    i<-grep(paste0("^",genename,"$"),rownames(seu_obj),value=TRUE,ignore.case = FALSE)
    if(length(i)==0){
      cat("\33[31m",genename,"not match:",grep(paste0("^",genename,"$"),rownames(seu_obj),value=TRUE,ignore.case = TRUE,perl = TRUE),"\33[0m\n")
    }
    else print(genename)
  }
}

png("classfy.png", width = 2000, height = 500)
Seurat::DotPlot(seu_obj, features =unlist(markers), dot.scale = 10) + Seurat::RotatedAxis()
dev.off()

celltype <- c("Oligodendrocytes",#1
              "Immue Cells",#2
              "Neuron",#3
              "Vascular cells",#4
              "Oligodendrocytes",#5
              "Astrocyte",#6
              "Neuron",#7
              "Neuron",#8
              "Neuron",#9
              "Vascular cells",#10
              "Neuron",#11
              "Immue cells",#12
              "Neuron",#13
              "Unclassfied",#14
              "Neuron" #15
)

Seurat::Idents(seu_obj) <- seu_obj@meta.data$seurat_clusters
names(celltype) <- levels(seu_obj)
seu_obj <- Seurat::RenameIdents(seu_obj,celltype)
seu_obj@meta.data$celltype <- Seurat::Idents(seu_obj)
