# sonia cebrian camison
# 4/agosto/2025
# estructura culex pipiens
# this script is to summarise enrichment results

#---------------
# libraries
#--------------
# for summary of enrichmnet
library(readr)
library(dplyr)
library(openxlsx)
# for extracting list of GO terms
library(httr)
library(jsonlite)

#------------------------------------------
# paste results of BP, MF and CC for each inversions
#-------------------------------
# 1. Define input files and their category labels
inversion <- "3p"  # change as needed
categories <- c("BP", "MF", "CC")
files <- paste0(inversion, "_enrich_", categories, "_significant.tsv")

# 2. Read and label each table
tables <- lapply(seq_along(files), function(i) {
  read_tsv(files[i], show_col_types = FALSE) %>%
    mutate(Category = categories[i]) %>%
    select(Category, everything())
})

# 3. Combine all into one data frame
merged_table <- bind_rows(tables)

# 4. Write to Excel
write.xlsx(merged_table, paste0(inversion, "_GO_enrichment.xlsx"), rowNames = FALSE)

#------------------------
# list of GO terms candidates
# ------------------------
# i have a list of candidates with many GO per gene. Is difficult to review by hand all 
#Go terms so i will extract them here

inversion <- "1q"  # change as needed

# 1. Read file and extract all GO terms
lines <- readLines(paste0(inversion,  "_candidates_GO.tsv"))
go_terms <- unique(unlist(regmatches(lines, gregexpr("GO:\\d{7}", lines))))
go_terms <- sort(go_terms)

# 2. Define function to query QuickGO
get_go_name <- function(go_id) {
  url <- paste0("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/", go_id)
  res <- GET(url, accept("application/json"))
  if (status_code(res) == 200) {
    content(res, "parsed")$results[[1]]$name
  } else {
    NA
  }
}

# 3. Query all GO terms
go_df <- data.frame(GO_term = go_terms, Name = sapply(go_terms, get_go_name))

# 4. Save to file

write.xlsx(go_df, paste0(inversion, "_candidates_GO_dictionary.xlsx"), rowNames = FALSE)

