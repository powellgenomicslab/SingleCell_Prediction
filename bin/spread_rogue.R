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

output <- file.path("results", paste0(date, output_dir_name))

if(!dir.exists(output)){
  dir.create(here(output))
}

# Read data ---------------------------------------------------------------

# Input

input    <- file.path("data", "2018-10-18_rogue_cells") # <------ Input directory
filename <- "counts.txt" # <------ Input file


# Read file

data <- fread(input = here(input, filename), header = FALSE, data.table = FALSE)
data <- spread(data, key = "V3", 2, fill = 0) %>% column_to_rownames("V1")
write_delim(data, path = here(output, "rogue_matrix.txt"), delim = "\t")



# Session info ------------------------------------------------------------

options(width = 70)
capture.output(devtools::session_info(), file = here(output, "session_info.txt"))