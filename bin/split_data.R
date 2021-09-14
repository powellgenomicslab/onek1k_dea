#   ____________________________________________________________________________
#   Script information                                                      ####

# title: Split count OneK1K data by individual and cell type
# author: Jose Alquicira Hernandez
# date: 2021-09-13

#   ____________________________________________________________________________
#   HPC details                                                             ####

# screen -S dea
# qrsh -N dea -l mem_requested=150G -q short.q
# conda activate r-4.1.0

#   ____________________________________________________________________________
#   Import libraries                                                        ####

library("dsLib")
library("data.table")
library("Seurat")
library("SeuratExtra")
library("SeuratDisk")

#   ____________________________________________________________________________
#   Set output                                                              ####

output <- set_output("2021-09-13", "split-individual_celltype")

#   ____________________________________________________________________________
#   Import data                                                             ####

inicio("Read data")
data <- LoadH5Seurat(here("..",
                          "onek1k_azimuth", 
                          "results", 
                          "2021-06-10_cell_annotation",
                          "onek1k_seurat.h5seurat"), 
                     assays = list(RNA = "counts"))
Idents(data) <- "predicted.celltype.l2"
fin()


covariates <- fread(here("data", "age_sex_info.tsv"))
setnames(covariates, "sampleid", "individual")

#   ____________________________________________________________________________
#   Create sample metadata                                                  ####

md <- as.data.table(data[[]][, c("individual", "pool")])
md <- unique(md)

stopifnot(all(md$individual %in% covariates$individual) & 
            all(covariates$individual %in% md$individual))

sample_md <- merge(md, covariates, by = "individual")

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

n <- as.data.table(data[[]][, c("individual", "predicted.celltype.l2")])
n <- n[,.N, .(individual, predicted.celltype.l2)]
setnames(n, "predicted.celltype.l2", "cell_type")

sample_md_celltype <- lapply(res, colnames)
sample_md_celltype <- lapply(sample_md_celltype, function(x){ 
  x <- as.data.table(x)
  setnames(x, "individual")
  })
sample_md_celltype <- lapply(sample_md_celltype, merge, sample_md, "individual")

sample_md_celltype <- mapply(function(x, cell_type){
  x[, cell_type := ..cell_type]
  x[, pair := paste(individual, ..cell_type, sep = "-")]
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
  colnames(x) <- paste(colnames(x), cell_type, sep = "-")
  x
}, res, names(res))

#   ____________________________________________________________________________
#   Export data                                                             ####

fwrite(sample_md, here(output, "sample_metadata.csv"))
saveRDS(sample_md_celltype, here(output, "sample-celltype_metadata.RDS"))
saveRDS(res, here(output, "expression_sample-celltype.RDS"))

#   ____________________________________________________________________________
#   Session info                                                            ####

print_session(here(output))

