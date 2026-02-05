library(Seurat)
library(future)
plan("multisession", workers = 1)
library(ggplot2)
options(future.seed = TRUE)
options(future.globals.maxSize = 40000 * 1024^2)
library(spacexr)
set.seed(123)
xenium.obj = LoadXenium("/rsrch3/home/genomic_med/ssatpati/Amit_Moran/Xenium/SCC_BCC_FF_IS_IC/20250116__221632__SCC_BCC_JHJG99_01162025/output-XETG00074__0024835__813_SCC__20250116__221656/", fov = "fov")
query.counts <- GetAssayData(xenium.obj, assay = "Xenium", layer = "counts")
coords <- GetTissueCoordinates(xenium.obj, which = "centroids")
rownames(coords) <- coords$cell
coords$cell <- NULL

allen.cortex.ref <- readRDS("/rsrch3/home/genomic_med/ssatpati/Amit_Moran/Xenium/res_0.2/res0.2-to-subsetcluster1/All_SCC_BCC_Xenium_final_res0.2_clust1-subclustered.rds")
head(rownames(allen.cortex.ref))

allen_cells_stripped <- sub("^[A-Za-z0-9_]+_[0-9]+_", "", Cells(allen.cortex.ref))
xenium.obj_sub <- subset(xenium.obj, cells = allen_cells_stripped)
xenium.obj_sub$sub <- allen.cortex.ref@meta.data$sub[match(Cells(xenium.obj), allen_cells_stripped)]
head(xenium.obj_sub@meta.data)
keep_cells <- rownames(xenium.obj_sub@meta.data[!is.na(xenium.obj_sub@meta.data$sub), ])
xenium.obj_sub2 <- subset(xenium.obj_sub, cells = keep_cells)
head(xenium.obj_sub2@meta.data)

print("first")
head(rownames(xenium.obj_sub2))
head(Cells(xenium.obj_sub2))


xenium.obj_sub2 <- BuildNicheAssay(object = xenium.obj_sub2, fov = "fov", group.by = "sub",niches.k = 15, neighbors.k = 30)

color_palette <- c("#442288", "#6CA2EA", "#FF00FF", "#FED23F", "#EB7D5B",
                   "#FF5733", "#7FFF00", "#33DBFF", "#900C3F", "#AF46B4",
                   "#32CD32", "#B5D33D", "#0047AB", "#7FFFD4", "#3DB5D3",
                   "#85660D","#F8A19F","#1C8356")

celltype.plot <- ImageDimPlot(xenium.obj_sub2, group.by = "sub", size = 1.5, cols = "polychrome", dark.background = F) + ggtitle("Cell type")
niche.plot <- ImageDimPlot(xenium.obj_sub2, group.by = "niches", size = 1.5, dark.background = F) + ggtitle("Niches") +
  scale_fill_manual(values = color_palette[1:18])

combined_plot <- celltype.plot | niche.plot

ggsave("niche_plot.pdf", plot = combined_plot, device = "pdf", width = 10, height = 6)

table(xenium.obj_sub2$sub, xenium.obj_sub2$niches)

table_data <- table(xenium.obj_sub2$sub, xenium.obj_sub2$niches)

table_matrix <- as.matrix(table_data)

# Write the data frame to a CSV file
write.csv(table_matrix, file = "niches.csv", row.names = TRUE)
unique(xenium.obj_sub2$niches)

table(is.na(xenium.obj_sub2$niches))

head(xenium.obj_sub2@meta.data)

print("third")
head(rownames(xenium.obj_sub2))
# Save te Seurat object (e.g., tmp) to an RDS file
saveRDS(xenium.obj_sub2, file = "niche.rds")





