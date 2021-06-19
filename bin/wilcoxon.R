#   ____________________________________________________________________________
#   Script information                                                      ####

# title:
# author: Jose Alquicira Hernandez
# date: 2021-06-10
# description: None

#   ____________________________________________________________________________
#   HPC details                                                             ####

# screen -S dea
# qrsh -N dea -l mem_requested=150G -pe smp 5 -q short.q
# conda activate r-4.1.0

#   ____________________________________________________________________________
#   Import libraries                                                        ####

# Primary
library("data.table")
library("dsLib")

# Secondary
library("Seurat")
library("SeuratDisk")
library("future")


#   ____________________________________________________________________________
#   Define future                                                           ####

options(future.globals.maxSize = 150 * 1024^3)
plan(multicore(workers = 5))


#   ____________________________________________________________________________
#   Set output                                                              ####

output <- set_output("2021-06-10", "wilcoxon")

#   ____________________________________________________________________________
#   Import data                                                             ####

inicio("Read data")
data <- LoadH5Seurat(here("..",
                  "onek1k_azimuth", 
                  "results", 
                  "2021-06-10_cell_annotation",
                  "onek1k_seurat.h5seurat"), 
             assays = list(SCT = "data"), 
             reductions = "sumap")

Idents(data) <- "predicted.celltype.l2"
fin()

#   ____________________________________________________________________________
#   Find DEG using wilcoxon test                                            ####

inicio("Find markers")
markers_wilcox <- FindAllMarkers(data, logfc.threshold = 0.5, only.pos = TRUE)
fin()


markers_wilcox <- as.data.table(markers_wilcox)

#   ____________________________________________________________________________
#   Export data                                                             ####

fwrite(markers_wilcox, here(output, "markers_wilcoxon_logfc-0.5.csv"))

#   ____________________________________________________________________________
#   Session info                                                            ####

print_session(here(output))

