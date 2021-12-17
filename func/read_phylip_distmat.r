#' Read a phylip formatted distance matrix
#'
#' From https://rdrr.io/github/SWittouck/tidyorthogroups/src/R/readers.R
#'
#' This function reads a phylip formatted distance matrix, as generated by e.g.
#' the EMBOSS tool distmat.
#' 
#' The number of lines to skip can be different for different phylip distance
#' matrix files.
#' 
#' For each row, the value of sequence_1 will be before the value of sequence_2
#' in sorting order.
#'
#' @param path Path to a phylip formatted distance matrix file
#' @param skip Number of lines before the actual distance values start
#' @param include_diagonal Should the diagonal of the distance matrix be included? 
#' 
#' @return A tibble with the variables sequence_1, sequence_2 and distance
#' 
#' @export
read_phylip_distmat <- function(path, skip = 8, include_diagonal = T) {
  
  distances_raw <- 
    readr::read_tsv(path, col_names = F, skip = skip) %>%
    separate(ncol(.), into = c("name", "number"), sep = " ")
  
  names <- distances_raw$name
  
  n <- length(names)
  
  distances <- 
    distances_raw %>%
    select_if(~ ! all(is.na(.))) %>%
    select(1:(!! n)) %>%
    `names<-`(names) %>%
    mutate(sequence_1 = !! names) %>%
    gather(key = "sequence_2", value = "distance", - sequence_1, na.rm = T) %>%
    mutate_at("distance", as.double)
  
  distances_1 <- filter(distances, sequence_1 >= sequence_2)
  
  distances_2 <- 
    distances %>%
    filter(sequence_1 < sequence_2) %>%
    mutate(sequence_1_temp = sequence_2, sequence_2_temp = sequence_1) %>%
    select(sequence_1 = sequence_1_temp, sequence_2 = sequence_2_temp, distance)
  
  distances <- bind_rows(distances_1, distances_2)
  
  if (! include_diagonal) {
    distances <- filter(distances, sequence_1 != sequence_2)
  }
  
  distances
  
}
