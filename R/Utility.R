## Some Utility functions to facilitate data processing
## Copyright (C) 2018-2020 AHM Mahfuzur Rahman (rahma118@umn.edu)
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.


#' Read the interaction/fitness data given a filename (a usual tab delimited format or R style format)
#'
#' @param filename.input: Input data file; genes on the row and samples on the column. First column provides the gene names, and the rest provide numeric data (interaction, fitness, etc.).
#' @param delim: Delimiter used to separate columns (tab is the default)
#' @param filename.genes: A file with a set of genes (exclude anything else).
#' @param filename.queries: A file with a set of samples (exclude anything else).
#'
#' @return Interaction data (numeric) with genes as rownames() and queries as colnames()
#' Ensures that the interaction data is sorted alphabetically (by gene names)
#' @examples
#' file.int <- 'Ceres_score_19Q2_with_depmap_id.txt'
#' data.interaction <- GetInteractionData(file.int)
#' @export
#' 

GetInteractionData <- function(filename.input, delim = '\t', filename.genes = NULL, filename.queries = NULL){
  
  data.interaction <- read.table(filename.input, header=TRUE, sep=delim)
  names.queries <- colnames(data.interaction)
  
  # Option1: If the input file is a usual format (with column names for each columns)
  if (class(data.interaction[,1]) != 'numeric'){
    names.genes <- data.interaction[,1]
    data.interaction[,1] <- NULL # data.interaction <- data.interaction[,-c(1)]
    rownames(data.interaction) <- names.genes
  }
  # Option2: If the input file is more of an R format (no column names for first column)
  else {
    names.genes <- row.names(data.interaction)
  }
  
  # sort data by gene names (if not sorted already)
  if (is.unsorted(names.genes)){
    ind <- order(names.genes)
    names.genes <- names.genes[ind]
    data.interaction <- data.interaction[ind,]
  }
  
  if (!is.null(filename.genes)) {
    # If we only need a subset of names.genes
    names.genes.small.list <- read.table(filename.genes, header = FALSE, row.names = NULL)
    
    tmp.data.genes <- intersect(names.genes.small.list[,1], names.genes)
    data.interaction <- data.interaction[tmp.data.genes,]
  }
  
  if (!is.null(filename.queries)) {
    # If we only need a subset of queries
    names.queries.small.list <- read.table(filename.queries, header = FALSE, row.names = NULL)
    queries.tmp <- intersect(names.queries.small.list[,1], names.queries)
    data.interaction <- data.interaction[,queries.tmp]
  }
  
  return(data.interaction)
}


#' Reads the Co-annotation standard data
#'
#' @param filename Name of the file where the co-annotation data is. Expected: A file with three/four columns: gene1 (string) gene2 (string) co-annotation (0/1) source (IDs).
#' @param delim: Delimiter used to separate columns (tab is the default)
#'
#' @return co.annotation The co-annotation data
#' @export

GetCoAnnotationData <- function(filename, delim = '\t'){
  
  co.annotation <- read.table(filename, header=FALSE, sep=delim, stringsAsFactors = FALSE)
  colnames(co.annotation) <- c('gene1', 'gene2', 'is_annotated', 'ID')
  
  if (is.unsorted(co.annotation$gene1)){
    ind <- order(co.annotation$gene1)
    co.annotation <- co.annotation[ind,]
  }
  
  # If there are no explicit returns from a function, 
  # the value of the last evaluated expression is returned automatically in R.
  co.annotation
}



#' Read the DepMap data and convert it to a gene by screen format
#'
#' @param filename Name/location of the input depmap file.
#' @param delim: Delimiter used to separate columns (',' is the default: usually used in depmap scores)
#' 
#' @return a data.frame in the format of genes (row) * screens (columns)
#' 
#' @export
#' 
ReformatDepMapData <- function(filename, delim = ','){

  # filename1 <- 'gene_effect.csv'
  data <- read.table(filename, header = F, sep = delim, stringsAsFactors = F)
  
  # Get the genes
  genes_tmp <- as.character(data[1,])
  genes_tmp <- genes_tmp[-1]
  # https://stackoverflow.com/questions/31235165/sapply-with-strsplit-in-r
  genes <- sapply(strsplit(genes_tmp, ' '), '[', 1) # Now we need to remove the trailing parts after the gene name
  
  # Way 1: Need a bit modification
  # screens <- data[,1]
  # screens <- screens[-1]
  # data.interaction.final <- t(as.numeric(as.character(data[2:dim(data)[1], 2:dim(data)[2]])))
  
  # Way 2: Read the data again (it will be numeric only now)
  data.interaction <- read.table(filename, header = T, sep = ',', stringsAsFactors = F, row.names = 1)
  data.interaction.final <- t(data.interaction)
  row.names(data.interaction.final) <- genes  
  
  return(as.data.frame(data.interaction.final))
}



#' Calculate the unique elements and their ranges
#'
#' @param data a sorted vector
#'
#' @return group.indices range of indices of unique elements (with rownames as the unique elements themselves)
#' 
GroupUniqueElements <- function(data){
  
  # gives the first occurrence of the data -> logical index
  logical.ind.unique <- !duplicated(data)
  elements.unique <- data[logical.ind.unique] # unique(data) 
  unique.index <- which(!duplicated(data)) # will give the indices of the first unique element
  # system.time({unique(data)})
  # system.time({!duplicated(data)})
  
  group.indices <- matrix(0, length(unique.index), 2)
  
  # If there is only one element
  if (length(unique.index) == 1){
    group.indices[1,1] <- 1 
    group.indices[1,2] <- length(data)
  }
  
  for (i in 1 : length(unique.index) - 1){
    group.indices[i,1] <- unique.index[i]
    group.indices[i,2] <- unique.index[i+1]-1
  }
  
  # For the last element
  i <- i + 1
  group.indices[i,1] <- unique.index[i] 
  group.indices[i,2] <- length(data)
  
  # Add row and column names
  colnames(group.indices) <- c('start', 'end')
  row.names(group.indices) <- elements.unique
  
  # return (list(elements.unique = elements.unique, group.indices = group.indices))
  return (group.indices)
}
