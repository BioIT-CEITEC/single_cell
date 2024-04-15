library(Seurat)

seurat_object <- function(args){
  current_dir <- getwd()
  input_dir <- args[1]
  output_dir <- args[2]
  sample_name <- args[3]
  plot_dir <- file.path(output_dir, "Plots")
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  sce <- Seurat::ReadSTARsolo(input_dir)
  setwd(output_dir)
  seurat_starsolo_CR <- Seurat::CreateSeuratObject(sce,project = sample_name,min.cells = 2, min.features = 200)
  saveRDS(seurat_starsolo_CR, file = paste0(sample_name, ".rds"))
  
  seurat_starsolo_CR <- Seurat::NormalizeData(seurat_starsolo_CR)
  print("Normalization done")
  seurat_starsolo_CR <- Seurat::FindVariableFeatures(seurat_starsolo_CR, selection.method = "vst", nfeatures = 2000)
  print("Find variable features done")
  seurat_starsolo_CR <- Seurat::ScaleData(seurat_starsolo_CR, features = rownames(seurat_starsolo_CR)) # means all genes
  print("Scale data done")
  seurat_starsolo_CR <- Seurat::RunPCA(seurat_starsolo_CR, features = VariableFeatures(object = seurat_starsolo_CR))
  print("PCA done")
  seurat_starsolo_CR <- Seurat::RunTSNE(seurat_starsolo_CR)
  print("TSNE done")
  seurat_starsolo_CR <- Seurat::RunUMAP(seurat_starsolo_CR, dims = 1:20)
  print("UMAP done")

  #newer definition of cluster by resolution - unknown number of cluster on the end
  # seurat_starsolo_CR <- FindNeighbors(seurat_starsolo_CR, dims = 1:20)
  # seurat_starsolo_CR <- FindClusters(seurat_starsolo_CR, resolution = c(.05,.1,0.2,.3,0.5,1))
  #visualisation based on newer set up
  # DimPlot(seurat_starsolo_CR,reduction = "umap", group.by = "RNA_snn_res.0.1", label = TRUE)
  # DimPlot(seurat_starsolo_CR,reduction = "tsne", group.by = "RNA_snn_res.0.1", label = TRUE)

  setwd("Plots")
  # LOOP for creating 1 - 12 cluster for data
  for (clust_num in (1:12)){
    #k means clustering in PCA reduction data
    kmeans_obj <- kmeans(Embeddings(seurat_starsolo_CR, reduction = "pca"), centers = clust_num)
    print("kmeans_obj produced")
    # add cluster to metadata
    seurat_starsolo_CR <- AddMetaData(seurat_starsolo_CR, metadata = kmeans_obj$cluster, col.name = "ClusterID")
    print("seurat_starsolo_CR produced")

    # UMAP---
    umap_plot<-plot(Seurat::DimPlot(seurat_starsolo_CR, reduction = "umap", group.by = "ClusterID", label = T)+ggplot2::ggtitle(""))
    print("umap_plot object created")
    png(filename =  paste0(sample_name,"_",clust_num, "_UMAP.png"))
    #ggplot2::ggsave(umap_plot,filename =  paste0(sample_name,"_",clust_num, "_UMAP.png"),create.dir = T )
    print(umap_plot)
    dev.off()
    print("umap_plot saved")

    #TSNE---
    tsne_plot<-plot(Seurat::DimPlot(seurat_starsolo_CR, reduction = "tsne", group.by = "ClusterID", label = T)+ggplot2::ggtitle(""))
    
    print("tsne_plot object created")
    png(filename = paste0(sample_name,"_",clust_num, "_TSNE.png"))
    print(tsne_plot)
    dev.off()
   # ggplot2::ggsave(tsne_plot,filename = paste0(sample_name,"_",clust_num, "_TSNE.png"),create.dir = T )
    print("tsne_plot saved")
  }


  setwd(current_dir)
}

#run as Rscript
args <- commandArgs(trailingOnly = TRUE)
seurat_object(args)