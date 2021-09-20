#   ____________________________________________________________________________
#   Script information                                                      ####

# title: Split OneK1K data by cell types
# author: Jose Alquicira Hernandez
# date: 2021-09-20

#   ____________________________________________________________________________
#   HPC details                                                             ####

# screen -S dea
# qrsh -N dea -l mem_requested=150G -q short.q
# conda activate r-4.1.0

#   ____________________________________________________________________________
#   Import libraries                                                        ####

library("dsLib")
library("purrr")
library("Seurat")
library("SeuratDisk")

#   ____________________________________________________________________________
#   Set output                                                              ####

output <- set_output("2021-09-20", "split-celltype")

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

#   ____________________________________________________________________________
#   Split object                                                            ####

inicio("Split data")
data <- SplitObject(data)
fin()

#   ____________________________________________________________________________
#   Export data                                                             ####

iwalk(data, ~ saveRDS(.x, file = here(output, paste(.y, ".RDS"))))

#   ____________________________________________________________________________
#   Session info                                                            ####

print_session(here(output))

