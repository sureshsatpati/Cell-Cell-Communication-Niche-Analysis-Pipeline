library(Seurat)
library(future)
plan("multisession", workers = 10)
library(ggplot2)
options(future.seed = TRUE)
options(future.globals.maxSize = 40000 * 1024^2)
library(spacexr)

xenium.obj = LoadXenium("/rsrch3/home/genomic_med/ssatpati/Amit_Moran/Xenium/SCC_BCC_FF_IS_IC/20250123__213134__SCC_BCC_JHJG99_01232025/output-XETG00074__0025112__11_SCC__20250123__213152/", fov = "fov")
query.counts <- GetAssayData(xenium.obj, assay = "Xenium", layer = "counts")
coords <- GetTissueCoordinates(xenium.obj, which = "centroids")
rownames(coords) <- coords$cell
coords$cell <- NULL
query <- SpatialRNA(coords, query.counts, colSums(query.counts))

allen.cortex.ref <- readRDS("/rsrch3/home/genomic_med/ssatpati/Amit_Moran/Xenium/res_0.2/All_SCC_BCC_Xenium_final.rds")
Idents(allen.cortex.ref) <- "seurat_clusters"
counts <- GetAssayData(allen.cortex.ref, assay = "Xenium", layer = "counts")
cluster <- as.factor(allen.cortex.ref$seurat_clusters)
names(cluster) <- colnames(allen.cortex.ref)
nUMI <- allen.cortex.ref$nCount_Xenium
names(nUMI) <- colnames(allen.cortex.ref)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
reference <- Reference(counts, cluster, nUMI)

allen_cells_stripped <- sub("^[A-Za-z0-9_]+_[0-9]+_", "", Cells(allen.cortex.ref))
head(allen_cells_stripped)
head(Cells(xenium.obj))
allen.cortex.ref@meta.data$cleaned_cell_id <- allen_cells_stripped
xenium.obj$seurat_clusters <- allen.cortex.ref@meta.data$seurat_clusters[match(Cells(xenium.obj), allen_cells_stripped)]
keep.cells <- Cells(xenium.obj)[!is.na(xenium.obj$seurat_clusters)]
xenium.obj_sub <- subset(xenium.obj, cells = keep.cells)
head(xenium.obj_sub)
xenium.obj_sub <- BuildNicheAssay(object = xenium.obj_sub, fov = "fov", group.by = "seurat_clusters",niches.k = 18, neighbors.k = 30)
celltype.plot <- ImageDimPlot(xenium.obj_sub, group.by = "seurat_clusters", size = 1.5, cols = "polychrome", dark.background = F) + ggtitle("Cell type")

# Check the number of levels in 'niches'
#n_levels <- length(levels(xenium.obj_sub$niches))

# Create a color palette with the appropriate number of colors
#colors <- c("#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B")[1:n_levels]

color_palette <- c("#442288", "#6CA2EA", "#FF00FF", "#FED23F", "#EB7D5B",
                   "#FF5733", "#7FFF00", "#33DBFF", "#900C3F", "#AF46B4",
                   "#32CD32", "#B5D33D", "#0047AB", "#7FFFD4", "#3DB5D3",
                   "#85660D","#F8A19F","#1C8356")

# Then apply the palette
#niche.plot <- ImageDimPlot(xenium.obj_sub, group.by = "niches", size = 1.5, dark.background = F) + 
#  ggtitle("Niches") + 
#  scale_fill_manual(values = c("#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B"))

niche.plot <- ImageDimPlot(xenium.obj_sub, group.by = "niches", size = 1.5, dark.background = F) + ggtitle("Niches") +
  scale_fill_manual(values = color_palette[1:18])



combined_plot <- celltype.plot | niche.plot

ggsave("combined_plot.pdf", plot = combined_plot, device = "pdf", width = 10, height = 6)

table(xenium.obj_sub$seurat_clusters, xenium.obj_sub$niches)

table_data <- table(xenium.obj_sub$seurat_clusters, xenium.obj_sub$niches)

table_matrix <- as.matrix(table_data)
#write.csv(master_niche_mat,"master_niche_mat_scc_ic_is.csv")

#df <- data.frame(seurat_clusters = xenium.obj_sub$seurat_clusters, niches = xenium.obj_sub$niches)

# Write the data frame to a CSV file
write.csv(table_matrix, file = "xenium_celltype_niches.csv", row.names = FALSE)

