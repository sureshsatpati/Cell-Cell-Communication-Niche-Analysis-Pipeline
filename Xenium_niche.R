#/rsrch4/home/genomic_med/METI/Suresh/Xenium/LPS
library(dplyr)
library(Seurat)
library(future)
plan("multisession", workers =20)
library(ggplot2)
options(future.seed = TRUE)
options(future.globals.maxSize = 40000 * 1024^2)

#list.object <- list(
#  xenium.obj.366 = CreateSeuratObject(counts = Read10X_h5("/rsrch3/home/genomic_med/ssatpati/Amit_Moran/Xenium/SCC_BCC_FF_IS_IC/20250109__220903__SCC_BCC_JHJG99_01092025/output-XETG00074__0036847__336_SCC__20250109__220931/cell_feature_matrix.h5")$`Gene Expression`,assay = 'Spatial',min.cells = 10,min.features = 3),
#  xenium.obj.23 = CreateSeuratObject(counts = Read10X_h5("/rsrch3/home/genomic_med/ssatpati/Amit_Moran/Xenium/SCC_BCC_FF_IS_IC/20250109__220903__SCC_BCC_JHJG99_01092025/output-XETG00074__0036847__23_SCC__20250109__220931/cell_feature_matrix.h5")$`Gene Expression`,assay = 'Spatial',min.cells = 10,min.features = 3),
#  xenium.obj.218 = CreateSeuratObject(counts = Read10X_h5("/rsrch3/home/genomic_med/ssatpati/Amit_Moran/Xenium/SCC_BCC_FF_IS_IC/20250109__220903__SCC_BCC_JHJG99_01092025/output-XETG00074__0036847__218_SCC__20250109__220931/cell_feature_matrix.h5")$`Gene Expression`,assay = 'Spatial',min.cells = 10,min.features = 3),
#  xenium.obj.335 = CreateSeuratObject(counts = Read10X_h5("/rsrch3/home/genomic_med/ssatpati/Amit_Moran/Xenium/SCC_BCC_FF_IS_IC/20250109__220903__SCC_BCC_JHJG99_01092025/output-XETG00074__0036927__35_BCC__20250109__220931/cell_feature_matrix.h5")$`Gene Expression`,assay = 'Spatial',min.cells = 10,min.features = 3),
#  xenium.obj.654 = CreateSeuratObject(counts = Read10X_h5("/rsrch3/home/genomic_med/ssatpati/Amit_Moran/Xenium/SCC_BCC_FF_IS_IC/20250109__220903__SCC_BCC_JHJG99_01092025/output-XETG00074__0036927__654_BCC__20250109__220931/cell_feature_matrix.h5")$`Gene Expression`,assay = 'Spatial',min.cells = 10,min.features = 3),
#  xenium.obj.564 = CreateSeuratObject(counts = Read10X_h5("/rsrch3/home/genomic_med/ssatpati/Amit_Moran/Xenium/SCC_BCC_FF_IS_IC/20250116__221632__SCC_BCC_JHJG99_01162025/output-XETG00074__0024832__564_BCC__20250116__221656/cell_feature_matrix.h5")$`Gene Expression`,assay = 'Spatial',min.cells = 10,min.features = 3),
#  xenium.obj.723 = CreateSeuratObject(counts = Read10X_h5("/rsrch3/home/genomic_med/ssatpati/Amit_Moran/Xenium/SCC_BCC_FF_IS_IC/20250116__221632__SCC_BCC_JHJG99_01162025/output-XETG00074__0024832__723_SCC__20250116__221656/cell_feature_matrix.h5")$`Gene Expression`,assay = 'Spatial',min.cells = 10,min.features = 3),
#  xenium.obj.738 = CreateSeuratObject(counts = Read10X_h5("/rsrch3/home/genomic_med/ssatpati/Amit_Moran/Xenium/SCC_BCC_FF_IS_IC/20250116__221632__SCC_BCC_JHJG99_01162025/output-XETG00074__0024835__738_SCC__20250116__221656/cell_feature_matrix.h5")$`Gene Expression`,assay = 'Spatial',min.cells = 10,min.features = 3),
#  xenium.obj.813 = CreateSeuratObject(counts = Read10X_h5("/rsrch3/home/genomic_med/ssatpati/Amit_Moran/Xenium/SCC_BCC_FF_IS_IC/20250116__221632__SCC_BCC_JHJG99_01162025/output-XETG00074__0024835__813_SCC__20250116__221656/cell_feature_matrix.h5")$`Gene Expression`,assay = 'Spatial',min.cells = 10,min.features = 3),
#  xenium.obj.4 = CreateSeuratObject(counts = Read10X_h5("/rsrch3/home/genomic_med/ssatpati/Amit_Moran/Xenium/SCC_BCC_FF_IS_IC/20250123__213134__SCC_BCC_JHJG99_01232025/output-XETG00074__0025119__4_SCC__20250123__213151/cell_feature_matrix.h5")$`Gene Expression`,assay = 'Spatial',min.cells = 10,min.features = 3),
#  xenium.obj.6 = CreateSeuratObject(counts = Read10X_h5("/rsrch3/home/genomic_med/ssatpati/Amit_Moran/Xenium/SCC_BCC_FF_IS_IC/20250123__213134__SCC_BCC_JHJG99_01232025/output-XETG00074__0025119__6_BCC__20250123__213152/cell_feature_matrix.h5")$`Gene Expression`,assay = 'Spatial',min.cells = 10,min.features = 3),
 # xenium.obj.11 = CreateSeuratObject(counts = Read10X_h5("/rsrch3/home/genomic_med/ssatpati/Amit_Moran/Xenium/SCC_BCC_FF_IS_IC/20250123__213134__SCC_BCC_JHJG99_01232025/output-XETG00074__0025112__11_SCC__20250123__213152/cell_feature_matrix.h5")$`Gene Expression`,assay = 'Spatial',min.cells = 10,min.features = 3)
#)

#list.fov <- list(
xenium.fov.11 = LoadXenium("/rsrch3/home/genomic_med/ssatpati/Amit_Moran/Xenium/SCC_BCC_FF_IS_IC/20250123__213134__SCC_BCC_JHJG99_01232025/output-XETG00074__0025112__11_SCC__20250123__213152", fov = "fov")
#  xenium.fov.6 = LoadXenium("/rsrch3/home/genomic_med/ssatpati/Amit_Moran/Xenium/SCC_BCC_FF_IS_IC/20250123__213134__SCC_BCC_JHJG99_01232025/output-XETG00074__0025119__6_BCC__20250123__213152", fov = "fov"),
#  xenium.fov.4 = LoadXenium("/rsrch3/home/genomic_med/ssatpati/Amit_Moran/Xenium/SCC_BCC_FF_IS_IC/20250123__213134__SCC_BCC_JHJG99_01232025/output-XETG00074__0025119__4_SCC__20250123__213151", fov = "fov"),
#  xenium.fov.813 = LoadXenium("/rsrch3/home/genomic_med/ssatpati/Amit_Moran/Xenium/SCC_BCC_FF_IS_IC/20250116__221632__SCC_BCC_JHJG99_01162025/output-XETG00074__0024835__813_SCC__20250116__221656", fov = "fov"),
#  xenium.fov.738 = LoadXenium("/rsrch3/home/genomic_med/ssatpati/Amit_Moran/Xenium/SCC_BCC_FF_IS_IC/20250116__221632__SCC_BCC_JHJG99_01162025/output-XETG00074__0024835__738_SCC__20250116__221656", fov = "fov"),
#  xenium.fov.723 = LoadXenium("/rsrch3/home/genomic_med/ssatpati/Amit_Moran/Xenium/SCC_BCC_FF_IS_IC/20250116__221632__SCC_BCC_JHJG99_01162025/output-XETG00074__0024832__723_SCC__20250116__221656", fov = "fov"),
#  xenium.fov.564 = LoadXenium("/rsrch3/home/genomic_med/ssatpati/Amit_Moran/Xenium/SCC_BCC_FF_IS_IC/20250116__221632__SCC_BCC_JHJG99_01162025/output-XETG00074__0024832__564_BCC__20250116__221656", fov = "fov"),
#  xenium.fov.654 = LoadXenium("/rsrch3/home/genomic_med/ssatpati/Amit_Moran/Xenium/SCC_BCC_FF_IS_IC/20250109__220903__SCC_BCC_JHJG99_01092025/output-XETG00074__0036927__654_BCC__20250109__220931", fov = "fov"),
#  xenium.fov.335 = LoadXenium("/rsrch3/home/genomic_med/ssatpati/Amit_Moran/Xenium/SCC_BCC_FF_IS_IC/20250109__220903__SCC_BCC_JHJG99_01092025/output-XETG00074__0036927__35_BCC__20250109__220931", fov = "fov"),
#  xenium.fov.218 = LoadXenium("/rsrch3/home/genomic_med/ssatpati/Amit_Moran/Xenium/SCC_BCC_FF_IS_IC/20250109__220903__SCC_BCC_JHJG99_01092025/output-XETG00074__0036847__218_SCC__20250109__220931", fov = "fov"),
#  xenium.fov.23 = LoadXenium("/rsrch3/home/genomic_med/ssatpati/Amit_Moran/Xenium/SCC_BCC_FF_IS_IC/20250109__220903__SCC_BCC_JHJG99_01092025/output-XETG00074__0036847__23_SCC__20250109__220931", fov = "fov"),
#  xenium.fov.336 = LoadXenium("/rsrch3/home/genomic_med/ssatpati/Amit_Moran/Xenium/SCC_BCC_FF_IS_IC/20250109__220903__SCC_BCC_JHJG99_01092025/output-XETG00074__0036847__336_SCC__20250109__220931", fov = "fov")
#)

#xenium.obj.11 <- subset(xenium.obj.11, subset = nCount_Xenium > 0)
#xenium.obj.6 <- subset(xenium.obj.6, subset = nCount_Xenium > 0)
#xenium.obj.4 <- subset(xenium.obj.4, subset = nCount_Xenium > 0)
#xenium.obj.813 <- subset(xenium.obj.813, subset = nCount_Xenium > 0)
#xenium.obj.738 <- subset(xenium.obj.738, subset = nCount_Xenium > 0)
#xenium.obj.723 <- subset(xenium.obj.723, subset = nCount_Xenium > 0)
#xenium.obj.564 <- subset(xenium.obj.564, subset = nCount_Xenium > 0)
#xenium.obj.654 <- subset(xenium.obj.654, subset = nCount_Xenium > 0)
#xenium.obj.335 <- subset(xenium.obj.335, subset = nCount_Xenium > 0)
#xenium.obj.218 <- subset(xenium.obj.218, subset = nCount_Xenium > 0)
#xenium.obj.23 <- subset(xenium.obj.23, subset = nCount_Xenium > 0)
#xenium.obj.336 <- subset(xenium.obj.336, subset = nCount_Xenium > 0)

BuildNicheAssay.multiple_FOVs.MiniBatchKmeans <- function(
list.object,
list.fov,
group.by,
assay = "niche",
cluster.name = "niches",
neighbors.k = 20,
niches.k.range = 2:30 ,
batch_size = 20,
num_init = 20
) {
# check for fov in sample set
# remove if not found in object
remove = NULL # init list of indices to remove
for ( i in seq_along(list.object) ){ # message(i)
# get object and fov for each object
object = list.object[[i]]
fov = list.fov[[i]]

    if( fov %!in%  names(object@images) ){
        warning( "fov is not found in the i-th object.  Removing the object from the list.object and list.fov.  i =", i)
       remove = c(remove, i)
    }
}
for (i in rev(remove) ){
    list.object[[i]] = NULL
    list.fov[[i]] = NULL
}


for ( i in seq_along(list.object) ){ message(i)
    # get object and fov for each object
    object = list.object[[i]]
    fov = list.fov[[i]]
    
    # initialize an empty cells x groups binary matrix
    cells <- Cells( object[[fov]] )
    group.labels <- unlist(object[[group.by]][cells, ] )
    groups <- sort( unique(group.labels) )
    cell.type.mtx <- matrix(
        "data" = 0
        , "nrow" = length(cells)
        , "ncol" = length(groups)
    )
    rownames(cell.type.mtx) <- cells
    colnames(cell.type.mtx) <- groups

    # populate the binary matrix 
    cells.idx <- seq_along(cells)
    group.idx <- match(group.labels, groups)
    cell.type.mtx[cbind(cells.idx, group.idx)] <- 1
    
    # find neighbors based on tissue position
    coords <- Seurat::GetTissueCoordinates( object[[fov]], "which" = "centroids" )
    rownames(coords) <- coords[["cell"]]
    coords <- as.matrix(coords[ , c("x", "y")])
    neighbors <- Seurat::FindNeighbors(
        "object" = coords
        , "k.param" = neighbors.k # Defines k for the k-nearest neighbor algorithm
        , "compute.SNN" = F
    )
    
    # create niche assay
    sum.mtx <- as.matrix( neighbors[["nn"]] %*% cell.type.mtx )
    niche.assay <- CreateAssayObject( "counts" = t(sum.mtx) )
    object[[assay]] <- niche.assay
    DefaultAssay(object) <- assay
    
    # scale data 
    object <- ScaleData(object)
    
    # return edited object to list
    list.object[[i]] = object
    
    
}

# get aggregate data for ClusterR::MiniBatchKmeans
# columns = features
# rows = cells
# cells = values
DAT = lapply( seq_along(list.object), function(i){
    t( list.object[[i]][[assay]]@scale.data )
}  )
DAT <- do.call("rbind", DAT)



res.clusters = data.frame(row.names = rownames(DAT))

for ( k in niches.k.range ){ message("k=", k)
    # new column name
    newCol = paste0("kmeans_", k)
    # get centroids
    km_mb = ClusterR::MiniBatchKmeans(
        "data" = DAT
        , "clusters" = k # the number of clusters
        , "batch_size" = batch_size # the size of the mini batches
        , "num_init" = num_init # number of times the algorithm will be run with different centroid seeds
        , "max_iters" = 100 # the maximum number of clustering iterations. 
        , "init_fraction" = 0.2 # percentage of data to use for the initialization centroids (applies if initializer is kmeans++ or optimal_init). Should be a float number between 0.0 and 1.0.
        , "initializer" = "kmeans++" # the method of initialization. One of, optimal_init, quantile_init, kmeans++ and random. See details for more information
        , "early_stop_iter" = 10 # continue that many iterations after calculation of the best within-cluster-sum-of-squared-error
        , "verbose" = F
        , "CENTROIDS" = NULL
        , "tol" = 1e-04
        , "tol_optimal_init" = 0.3
        , "seed" = 1
    )
    
    # use centroids to get clusters

    res.clusters[,newCol] = ClusterR::predict_MBatchKMeans( # This function takes the data and the output centroids and returns the clusters.
        "data" = DAT
        , "CENTROIDS" = km_mb$centroids
    )
    res.clusters[,newCol] = as.factor( res.clusters[,newCol] ) # change clusters to factors
    
}

# get clusters back onto the objects
colnames(res.clusters) = paste0(cluster.name,".", colnames(res.clusters))
for ( i in seq_along(list.object) ){ message(i)
    # get object and fov for each object
    object = list.object[[i]]
    
    # get clusters in correct cell row order into metadata of object
    object[[]] = res.clusters[rownames(object[[]]),]
    
    # return edited object to list
    list.object[[i]] = object
}


return(list.object)
}


