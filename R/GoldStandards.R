## Functions to generate Co-annotation Gold Standards for FLEX
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


#' Create a co-annotation standard for gene pairs from given gene memberships into entities (complex, pw etc.).
#' 
#' @param data_standard the functional standard as a data.frame ($ID, $Name, $Genes as columns)
#' @param overlap_length Number of common annotations to a complex/pathway to qualify as a co-annoation
#' @param subset_str Set this (a vector of strings) if interested in a subset of the data (i.e. 'KEGG' for pathway)
#' @param file_location The location of the file (path + filename) to store the co-annotation data. If NULL is given, we still return the co_annotation data.
#'
#' @return a data frame for co-annotations of gene pairs
#'  gene1: combined.true: true co-ann values (from the standard 0/1) for pairs
#'  gene2: combined.score: predicted scores (according to similarity) for pairs
#'  is_annotated: 0/1 value (Are the gene pairs co-annotated to the standard)
#'  ID: Source(s) of the co-annotation
#'  
#' @examples 
#'  
#' @export
#'  

MakeCoAnnotationFromGeneSymbols <- function(data_standard, overlap_length = 1, subset_str = NULL, file_location = NULL){
  
  if (!is.null(file_location)){ # If a path for the file is given
    if (file.exists(file_location)){
      message('File already exists .. Loading data from file ...')
      
      load(file_location)
      ifelse (grep('data_co_annotation', ls()), return (data_co_annotation), return(NULL))
      
    } else{
      message("File doesn't exist! Will create a new one!!!")
    }
  }
  
  ## Subset the input data (only keep the Names that include the string as provided in subset_str)
  ind <- c()
  if (!is.null(subset_str)){
    for (i in 1 : length(subset_str)){
      tmp_ind <- grep(subset_str[i], data_standard$Name, ignore.case = FALSE)
      ind <- union(ind, tmp_ind)
    }
    data_standard <- data_standard[ind,]
  }
  
  
  ## Make a gene to ID (complex/pathway etc.) map -------------------------
  # https://blog.dominodatalab.com/a-quick-benchmark-of-hashtable-implementations-in-r/
  gene_id_map <- new.env() # env is better than other hash implementations in R (pretty fast, faster than Matlab containers.Map)
  get_hash <- Vectorize(get, vectorize.args = "x") # https://www.r-bloggers.com/hash-me-if-you-can/
  
  print('Mapping Genes to the standard (complex for example) IDs ... ')
  
  for (i in 1 : length(data_standard$ID)){
    # gene_list = unlist(strsplit(toupper(data_standard$Genes[i]), ';'))
    gene_list = unlist(strsplit(data_standard$Genes[i], split = ';')) # There are genes like c4orf3 where it may make difference
    gene_list = gsub(' ', '', gene_list) # Replacing any spaces with nothing
    
    for (j in 1 : length(gene_list)){
      key <- gene_list[j]
      key <- gsub(' ', '', key) # Remove spaces in name (if any)
      
      if (key == '') { # If it's an empty character
        next
      }
      
      # If the key is absent insert, otherwise append (the complex ID)
      # exists('MAP3K20', envir = gene_id_map) retruns T/F
      ifelse(is.null(gene_id_map[[key]]), gene_id_map[[key]] <- data_standard$ID[i], 
             gene_id_map[[key]] <- c(gene_id_map[[key]], data_standard$ID[i]))
    }
  }
  
  all_genes_in_std <- ls(gene_id_map) # names(gene_id_map) work too, but ls() outputs are sorted alphabetically
  entity_list_of_genes <- get_hash(all_genes_in_std, gene_id_map) # A list
  
  
  ## Now make the co-annotation standard -------------------------
  candidate_genes <-all_genes_in_std
  candidate_IDs <- entity_list_of_genes # A list with candiate_genes as list names
  n <- length(candidate_genes)
  
  gene_first <- rep('', choose(n, 2))
  gene_second <- rep('', choose(n, 2))
  annotation_value <- rep(0, choose(n, 2))
  complex_id_association <- rep('', choose(n, 2))
  count <- 1
  
  print('Generating Co-annotation Gold Standard ... ')
  pb <- txtProgressBar(style = 3) # Progress bar
  
  for (i in 1: (n-1)){
    for (j in (i+1):n){
      gene_first[count] <- candidate_genes[i]
      gene_second[count] <- candidate_genes[j]
      
      tmp <- intersect(candidate_IDs[[i]], candidate_IDs[[j]]) # test case, i = 1272 (HDAC1), j = 1274 (HDAC2)
      
      if (length(tmp) > 0){#'  
        ifelse(length(tmp) > 1, complex_id_association[count] <- paste0(as.character(tmp), collapse = ';'), 
               complex_id_association[count] <- as.character(tmp))  
      }
      
      annotation_value[count] <- length(tmp)
      count <- count + 1
    }
    setTxtProgressBar(pb, count/length(gene_first)) # Progress bar update
  }
  
  close(pb)
  
  # Convert annotation_value to a zero/one (presence or absence of co-annotation)
  data_co_annotation <- data.frame(gene1 = gene_first, gene2 = gene_second, is_annotated = annotation_value, 
                                   ID = complex_id_association, 
                                   stringsAsFactors = FALSE) # stringAsFactors is super important here
  data_co_annotation$is_annotated[data_co_annotation$is_annotated < overlap_length] <- 0
  data_co_annotation$is_annotated[data_co_annotation$is_annotated >= overlap_length] <- 1

  # If the direcotry itself doesn't exist, we still return the output
  if (!is.null(file_location)){
    result = tryCatch(save(data_co_annotation, file = file_location), 
                      error=function(e) {
                        message(paste("Directory/File does not to exist / no permission", file_location))
                        # message("Here's the original error message:")
                        # message(e)
                        
                        # Choose a return value in case of error
                        return(data_co_annotation)
                      },
                      finally = {})
  }
  
  return(data_co_annotation)
}



## ====================== This part is Work in Progress (and thus commented out) ====================== ##
# Create a co-annotation standard for gene pairs from given gene memberships into Gene Ontology (GO) Biological Processes (BP).
# 
# @param data_standard the functional standard as a data.frame ($ID, $Name, $Genes as columns)
# @param overlap_length Number of common annotations to a complex/pathway to qualify as a co-annoation; for GO, we need a range. The default minimum is 10 and maximum is 300.
# @return a data frame for co-annotations of gene pairs
#  gene1: combined.true: true co-ann values (from the standard 0/1) for pairs
#  gene2: combined.score: predicted scores (according to similarity) for pairs
#  is_annotated: 0/1 value (Are the gene pairs co-annotated to the standard)
#  ID: Source(s) of the co-annotation
#  
#  @examples 
#  
# MakeCoAnnotationForGOBP <- function(data_standard, overlap_length = c(10, 300)){
# 
#   ## Which packages would be the best suited for this?
#   # topGO or GO.db or biomaRt?
#   
#   ## GO.db Package (issue: this package doesn't give us associations)
#   # The purpose of this package is to provide detailed
#   # information about the latest version of the Gene Ontologies. This package is updated biannually.
#   # https://bioconductor.org/packages/release/data/annotation/manuals/GO.db/man/GO.db.pdf
#   install.packages('GO.db') 
#   
#   # Use select() # ?select
#   
#   # GOBPANCESTOR Annotation of GO Identifiers to their Biological Process Ancestors
#   # The format is an R object mapping the GO BP terms to all ancestor terms, where an ancestor term
#   # is a more general GO term that precedes the given GO term in the DAG (in other words, the parents,
#   #                                                                      and all their parents, etc.)
#   # xx <- as.list(GOBPANCESTOR)
#   # GOBPCHILDREN, GOBPOFFSPRING, GOBPPARENTS (and similar for CC, and MFs)
#   
#   
# }
