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



#' Download the GIANT functional network and process it to use for FLEX
#'
#' @param file_location The location save the input GIANT file. If it is not there, it will be downloaded.
#'
#' @return a data frame for co-annotations of gene pairs
#'  gene1: Name of the first gene in the gene pair
#'  gene2: Name of the second gene in the gene pairs
#'  is_annotated: 0/1 value (Are the gene pairs functionally related)
#' 
#' @importFrom org.Hs.eg.db org.Hs.egSYMBOL
#' @importFrom AnnotationDbi mappedkeys
#'
#' @export
#'  
#' @examples
#' file_loc <- '/project/chadm/Mahfuz/CRISPR/2_HAP1/GIN_Analysis/GIN_Paper_Analysis/3_Standard_functional_evaluation/1_PR_Analysis/Data/global_standard.txt'
#' 	data.ca <- MakeFuncNetFromGIANT(file_loc)
#' 
MakeFuncNetFromGIANT <- function(file_location = NULL){
  # file_location <- '/home/mahfuz/Desktop/CRISPR/FLEX/R/Tests/Package_Test/Vignette_test/global_top_edges.txt'
  # file_location <- '/project/chadm/Mahfuz/Database/Greene_et_al/Func_net/global_top.gz'
  
  if (!file.exists(file_location)){
    message('Downloading GIANT data ...')
    
    # The following is the actual link of the file (located in: https://hb.flatironinstitute.org/download/global)
    # The best way would be to download data using a RESTful API (https://hb.flatironinstitute.org)
    download.file(url='https://s3-us-west-2.amazonaws.com/humanbase/networks/global_top.gz', destfile = file_location, method='curl')
  }
  
  # Read the downloaded data (R can read a .gz file as a text file, if it indeed is a text file)
  data.GIANT.entrez <- read.table(file_location,  header = F, sep = '\t', quote = '', stringsAsFactors = F)
  
  ## @importFrom org.Hs.eg.db org.Hs.egSYMBOL
  
  ## --------------------------------------------------
  # Use an entrez to gene symbol mapping to convert this standard (in Entrez) to usable with FLEX (in Symbol)
  genes.entrezID <- sort(unique(union(data.GIANT.entrez[,1], data.GIANT.entrez[,2])))
  genes.entrez.str <- sapply(genes.entrezID, toString) # as.character(genes.entrezID)
  
  # https://stuff.mit.edu/afs/athena/software/r/current/arch/i386_linux26/lib/R/library/org.Hs.eg.db/html/org.Hs.egSYMBOL.html  
  x <- org.Hs.eg.db::org.Hs.egSYMBOL
  mapped_genes <- AnnotationDbi::mappedkeys(x)
  xx <- as.list(x[mapped_genes])
  
  # *** The rownames are entrez IDs and the values are gene symbols
  gene.entrez.to.symbol <- unlist(xx[genes.entrez.str])
  
  ## --------------------------------------------------
  #  Now convert the data from entrez ID format to gene format
  data.standard <- data.GIANT.entrez
  names(data.standard) <- c('gene1', 'gene2', 'is_annotated')
  
  # Sort the data by gene2 followed by gene1
  if (is.unsorted(data.standard$gene2)){
    ind <- order(data.standard$gene2)
    data.standard <- data.standard[ind,]
  }
  
  ind <- order(data.standard$gene1)
  data.standard <- data.standard[ind,]
  
  # Group the pairs by gene1
  gene.indices <- GroupUniqueElements(data.standard$gene1)
  
  # Now create a new data replacing the entrezIDs by gene symbols
  # The size should not exceed the original data (may reduce in absence of matching)
  gene1.symbol <- vector('character', dim(data.standard)[1])
  gene2.symbol <- vector('character', dim(data.standard)[1])
  is.annotated <- numeric(dim(data.standard)[1])
  
  curr.ind <- 1 # Track the progress
  
  for (i in row.names(gene.indices)){
    # i <- row.names(gene.indices)[1]
    
    # This entrez ID doesn't have an associated gene_symbol
    if (!(i %in% names(gene.entrez.to.symbol))){
      next
    }
    
    ## ---------- Get the symbol for first gene ----------
    std.first.gene <- gene.entrez.to.symbol[i]
    
    # Now get the corresponding second genes and their associations
    indices.i <- gene.indices[i,1] : gene.indices[i,2]
    sec.genes.entrez <- as.character(data.standard$gene2[indices.i])
    co.ann.values.cand <- data.standard$is_annotated[indices.i]
    
    ## ---------- Get the symbol for second genes ----------
    # Remove the genes that don't have an entrez to symbol mapping
    # sum(sec.genes.entrez %in% names(gene.entrez.to.symbol))
    sec.genes.entrez.cand <- sec.genes.entrez[sec.genes.entrez %in% names(gene.entrez.to.symbol)]
    std.second.genes <- gene.entrez.to.symbol[sec.genes.entrez.cand]
    
    ## ---------- Now get the corresponding annotation ----------
    common.genes <- intersect(sec.genes.entrez, sec.genes.entrez.cand)
    ind <- match(common.genes, sec.genes.entrez)
    co.ann.values <- co.ann.values.cand[ind]
    
    curr.size <- length(co.ann.values)
    ind.to.fill <- curr.ind: (curr.ind + curr.size - 1)
    
    gene1.symbol [ind.to.fill] <- rep(std.first.gene, length(ind.to.fill))
    gene2.symbol [ind.to.fill] <- std.second.genes
    is.annotated [ind.to.fill] <- co.ann.values
    
    curr.ind <- curr.ind + curr.size # Update
  }
  
  
  data.GIANT <- data.frame(gene1 = gene1.symbol, gene2 = gene2.symbol, is_annotated = is.annotated)
  data.GIANT <- data.GIANT[-(curr.ind : length(is.annotated)), ]
  
  # If there is any NA in the gene names (which shouldn't be the case)
  ind.na <- which(is.na(data.GIANT$gene2) | is.nan(data.GIANT$gene2))
  if (length(ind.na) > 0){
    data.GIANT <- data.GIANT[-ind.na, ]
  }
  
  return (data.GIANT)
}
