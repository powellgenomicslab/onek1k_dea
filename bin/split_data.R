#   ____________________________________________________________________________
#   Script information                                                      ####

# title: Split count OneK1K data by individual and cell type
# author: Jose Alquicira Hernandez
# date: 2021-09-13

#   ____________________________________________________________________________
#   HPC details                                                             ####

# screen -S pbulk
# qrsh -N pbulk -l mem_requested=150G -q short.q
# conda activate r-4.1.1

#   ____________________________________________________________________________
#   Import libraries                                                        ####

library("dsLib")
library("data.table")
library("Seurat")
library("SeuratDisk")
library("popdea")

#   ____________________________________________________________________________
#   Set output                                                              ####

output <- here("results", "2021-11-10_create_pseudobulk")
dir.create(output)

#   ____________________________________________________________________________
#   Import data                                                             ####

inicio("Read data")
data <- LoadH5Seurat(here("..", "onek1k_scd", "results", 
                          "2021-11-10_add_metadata", "onek1k.h5seurat"), 
                     assays = list(RNA = "counts"))
Idents(data) <- "predicted.celltype.l2"
fin()


# Get European only
data <- data[, !data$ethnic_outlier]

#   ____________________________________________________________________________
#   Create sample metadata                                                  ####

md <- as.data.table(data[[]][, c("individual", "pool", "age", "sex")])
md <- unique(md)

#   ____________________________________________________________________________
#   Split data by individual and cell type                                  ####

inicio()
groups <- split_matrix(data, by = c("individual", "predicted.celltype.l2"))
fin()

#   ____________________________________________________________________________
#   Aggregate UMIs by individual-cell type                                  ####

groups <- map_matrix(groups, Matrix::rowSums)

#   ____________________________________________________________________________
#   Tidy data                                                               ####

res <- reduce_matrix(groups)

#   ____________________________________________________________________________
#   Create metadata for each cell type                                      ####

# Get number of cells per individual across each cell ty[e]
n <- as.data.table(data[[]][, c("individual", "predicted.celltype.l2")])
n <- n[,.N, .(individual, predicted.celltype.l2)]
setnames(n, "predicted.celltype.l2", "cell_type")

# Retrieve individual ids present in each cell type
sample_md_celltype <- lapply(res, colnames)

# Create metadata template with individual ides for each cell type
sample_md_celltype <- lapply(sample_md_celltype, function(x){ 
  x <- as.data.table(x)
  setnames(x, "individual")
  })

# Merge individual-level metadata with cell-type individual templates
sample_md_celltype <- lapply(sample_md_celltype, merge, md, "individual")

# Add cell type label and number of cells per individual
sample_md_celltype <- mapply(function(x, cell_type){
  x[, cell_type := ..cell_type]
  x[, pair := paste(..cell_type, individual, sep = "-")]
  merge(x, n, by = c("individual", "cell_type"))
}, sample_md_celltype, names(sample_md_celltype), SIMPLIFY = FALSE)

# Validate order of cell types in sample metadata and expression data
stopifnot(all(names(sample_md_celltype) == names(res)))

# Order metadata according to expression data
sample_md_celltype <- mapply(function(metadata, exp){
  i <- colnames(exp)
  metadata[match(i, metadata$individual)]
}, sample_md_celltype, res, SIMPLIFY = FALSE)

# Verify order of individuals in metadata and expression data
stopifnot({all(mapply(function(metadata, exp){
  all(metadata$individual == colnames(exp))
}, sample_md_celltype, res))})


#   ____________________________________________________________________________
#   Rename individuals by cell type                                         ####

res <- mapply(function(x, cell_type){
  colnames(x) <- paste(cell_type, colnames(x), sep = ".")
  x
}, res, names(res))

#   ____________________________________________________________________________
#   Export data                                                             ####

fwrite(md, here(output, "sample_metadata.csv"))
saveRDS(sample_md_celltype, here(output, "sample-celltype_metadata.RDS"))
saveRDS(res, here(output, "expression_sample-celltype.RDS"))

#   ____________________________________________________________________________
#   Session info                                                            ####

print_session(here(output))

