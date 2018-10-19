# Script information ------------------------------------------------------

# title: Spread rogue cell data
# author: José Alquicira Hernández
# date: 2018-10-18
# description: None


# Import libraries --------------------------------------------------------

# Primary

library("tidyverse")
library("data.table")
library("here")

# Secondary


# Set output --------------------------------------------------------------


output_dir_name <- "rogue_matrix" # <------ Output directory

date <- format(Sys.Date(), format = "%Y-%m-%d_")
date <- "2018-10-18_"


output <- file.path("results", paste0(date, output_dir_name))

if(!dir.exists(output)){
  dir.create(here(output))
}

# Read data ---------------------------------------------------------------

# Input

input    <- file.path("data", "2018-10-18_rogue_cells") # <------ Input directory
filename <- "counts.txt" # <------ Input file


# Read file

data <- fread(input = here(input, filename), header = TRUE, data.table = FALSE)
data <- spread(data, key = "CELL_ID", 2, fill = 0) %>% column_to_rownames("gene_id")
saveRDS(data, here(output, "rogue_matrix.RDS"))


# Session info ------------------------------------------------------------

options(width = 70)
capture.output(devtools::session_info(), file = here(output, "session_info.txt"))