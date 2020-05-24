## Functions to associate pairwise scores with ground truth values.
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


#' Calculate the profile similarity between all possible pairs, pair up their true co-annotation value and source of the co-annoation.
#' 
#' @param data.standard: the functional standard in a format (Gene1 (chr), Gene2 (chr), 1/0, source of annotation (chr))
#' @param data.interaction: interaction scores (returned from GetInteractionData)
#'                   or pairwise.correlation of interaction data
#
#' @return a list of vectors
#'  true: combined.true: true co-ann values (from the standard 0/1) for pairs
#'  predicted: combined.score: predicted scores (according to similarity) for pairs
#'  source: origin of the pair (complex, pathway or GO ID)
#'  indices: exact location of the pair in the co-annoation standard (assumes sorted)
#'  
#' @examples 
#'  Complex.Pairs <- CalculatePredictionAndTrueOnLibraryProfiles (data.standard, pairwise.correlation)
#' @export
#' 

CalculatePredictionAndTrueOnLibraryProfiles <- function(data.standard, data.interaction){
  ## *** Check input data format
  true_classes_count <- (class(data.standard) == 'data.frame') + (class(data.interaction) == 'data.frame')
  if(true_classes_count < 2){
    stop ("Check the inputs again ... they are supposed to be data.frame's")
    # stopifnot(true_classes_count < 3)
    # print ("Check the inputs again ... they are supposed to be data.frame's")
    # return (NULL)
  }
  
  # 1. Pairwise data in the format of Gene1 Gene2 Similarity
  if (dim(data.interaction)[2] == 3) {
    
    # We want to check the data types for this format
    if ((typeof(data.interaction[1,1]) == 'character') & (typeof(data.interaction[1,2]) == 'character') & 
        (typeof(data.interaction[1,3]) == 'double')){
      print('Similarity data provided ... calling FromGenePairSimilarity ...')
      return (FromGenePairSimilarity(data.standard, data.interaction))
    } else{
      print('ill-formatted data .. quitting ...')
      return(NULL)
    }
  } 
  
  # 2. Pairwise correlation / similarity data provided as a square matrix
  else if ( (dim(data.interaction)[1] == dim(data.interaction)[2]) & 
       isSymmetric(as.matrix(data.interaction)) ) {
    # Check: 1. Square, 2. Symmetric
    print('Pairwise correlation/similarity matrix provided ...')
    
    # Sort the pairwise data by gene names (we need to do it on both row and column side)
    if (is.unsorted(rownames(data.interaction))){
      gene.names <- rownames(data.interaction) # Names are already associated to the rows and columns
      ind <- order(gene.names)
      data.interaction <- data.interaction[ind,ind]      
    }
    
    pairwise.correlation <- as.matrix(data.interaction)
    rm(data.interaction) # Not needed anymore (It's a 17K * 17K matrix here, so might be a bottleneck in low RAM pcs)
    gc()
    
    return (FromAllPairwiseCorrelation(data.standard, pairwise.correlation))
  } 
  
  # 3. Genes (row) * screens (column) matrix provided (GI or fitness, for example)
  else{ 
    # Calculate pairwise.correlation
    
    # Sort the interaction.data by gene names (library side)
    if (is.unsorted(rownames(data.interaction))){
      gene.names <- rownames(data.interaction)
      ind <- order(gene.names)
      data.interaction <- data.interaction[ind,]
      row.names(data.interaction) <- gene.names[ind]
    }
    
    print('calculating pairwise correlation ...')
    
    # Check if there are a lot of NaNs in data?
    percent.nan <- sum(is.nan(as.matrix(data.interaction))) / (dim(data.interaction)[1] * dim(data.interaction)[2])
    
    if (percent.nan > 0.1){ # If we have more than 20% NaNs in the data
      pairwise.correlation <- cor(t(data.interaction), use = 'pairwise.complete.obs', method = 'pearson') # Slow, and 
    }else{
      pairwise.correlation <- cor(t(data.interaction), use = 'complete.obs', method = 'pearson') # Fast  
    }
    
    ## *** TODO: Add inner product as the similarity for cases where we have a lot of NaNs or smaller profile size

    rm(data.interaction)
    return (FromAllPairwiseCorrelation(data.standard, pairwise.correlation))
  }
}


# Supporting function for CalculatePredictionAndTrueOnLibraryProfiles
FromAllPairwiseCorrelation <- function(data.standard, pairwise.correlation){
  
  # print('speed optimizations ...')
  length.corr <- dim(pairwise.correlation)[1]
  gene.symbol <- rownames(pairwise.correlation)
  
  # Sort the co annotation data (data.standard) if not already sorted
  if (is.unsorted(data.standard$gene1)){
    ind <- order(data.standard$gene1)
    data.standard <- data.standard[ind,]
  }
  
  out <- GroupUniqueElements(data.standard$gene1)
  gene.indices <- out
  unique.names.genes <- row.names(out)
  
  # Fill out similarity values and corresponding co-annotation (1/0) for each
  # common gene pairs between similarity and co-annotation
  combined.true <- numeric(dim(data.standard)[1]) # 0 - default/ initialization value
  combined.score <- numeric(dim(data.standard)[1])
  source <- vector('character', dim(data.standard)[1])
  indices.in.standard <- numeric(dim(data.standard)[1])
  
  ## Probable Speedup: Not use rownames directly, but find some numeric representation
  int.genes.in.std <- intersect(gene.symbol, unique.names.genes)
  ind.int.genes.in.std <- match(int.genes.in.std, gene.symbol)
  # identical(gene.symbol[ind.int.genes.in.std], int.genes.in.std) # sanity-check
  
  curr.ind <- 1 # Track the progress
  print('Mapping pairs against their correlation values ... ')
  
  # Progress bar
  pb <- txtProgressBar(style = 3) 
  pb_count <- 0
  
  for (i in ind.int.genes.in.std){
    
    # Query gene entrez ID
    # common.gene.ID <- gene.symbol[i] # gene symbol
    # curr.unique.gene <- intersect(common.gene.ID, unique.names.genes)
    
    curr.unique.gene <- gene.symbol[i]
    
    ## Get standard data for this names.genes (with others associated with it)
    ind.curr <- gene.indices[curr.unique.gene,1] : gene.indices[curr.unique.gene,2]
    std.second.genes <- data.standard[ind.curr, 2]
    co.ann.values <- data.standard[ind.curr, 3]
    names(co.ann.values) <- std.second.genes
    
    # *** Considering only Upper triangle (both datasets must be sorted)
    data.corr <- pairwise.correlation[i, (i+1) :length.corr] # parentheses is important
    sim.second.genes <- gene.symbol[(i+1) : length.corr]
    
    # system.time(sim.second.genes <- gene.symbol[(i+1) : length.corr]) # Faster
    # system.time(sim.second.genes <- gene.symbol[-(1:i)])
    
    ## Assign the true lables (1,0) from the functional standard
    common.genes <- intersect(sim.second.genes, std.second.genes)
    values.true <- co.ann.values[common.genes]
    values.predicted <- data.corr[common.genes]
    curr.size <- length(values.true)
    
    # print(paste(i, length(values.true), sep = ' '))
    
    if (curr.size > 0){ # R will produce an error otherwise
      combined.true[curr.ind: (curr.ind + curr.size - 1)]  <-  values.true
      combined.score[curr.ind: (curr.ind + curr.size - 1)] <-  values.predicted
      
      # If we have to store the IDs and indices too
      if (dim(data.standard)[2] == 4){
        co.ann.IDs <- data.standard[ind.curr, 4]
        names(co.ann.IDs) <- std.second.genes
        source[curr.ind: (curr.ind + curr.size - 1)] <- co.ann.IDs[common.genes]
        
        values.indices <- ind.curr
        names(values.indices) <- std.second.genes
        indices.in.standard[curr.ind: (curr.ind + curr.size - 1)] <- values.indices[common.genes]
      }
    }
    
    curr.ind <- curr.ind + curr.size # Update
    
    # Progress bar update
    pb_count <- pb_count + 1
    setTxtProgressBar(pb, pb_count/length(ind.int.genes.in.std))
  }
  close(pb)
  
  ## Remove the unnecessary part of the data
  combined.true <- combined.true[1: (curr.ind - 1)]
  combined.score <- combined.score[1: (curr.ind - 1)]
  
  ## Remove the na correlation values and corresponding interactions
  ind.na <- which(is.na(combined.score) | is.nan(combined.score))
  if (length(ind.na) > 0){
    combined.score <- combined.score[-ind.na]
    combined.true <- combined.true[-ind.na]  
  }
  
  # If data.standard is of dimension 4 (contains ID information)
  if (dim(data.standard)[2] == 4){
    source <- source[1: (curr.ind - 1)]
    indices.in.standard <- indices.in.standard[1: (curr.ind - 1)]
    
    if (length(ind.na) > 0){
      source <- source[-ind.na]
      indices.in.standard <- indices.in.standard[-ind.na]
    }
  }
  
  # Return outputs as a list
  if (dim(data.standard)[2] == 4){
    return(list(true = combined.true, predicted = combined.score, ID = source, index = indices.in.standard))  
  } else{
    return(list(true = combined.true, predicted = combined.score))
  }
}


# Supporting function for CalculatePredictionAndTrueOnLibraryProfiles
FromGenePairSimilarity <- function(data.standard, data.interaction){
  
  # Sort the co annotation data (data.standard)
  if (is.unsorted(data.standard$gene1)){
    ind <- order(data.standard$gene1)
    data.standard <- data.standard[ind,]
  }
  system.time(out <- GroupUniqueElements(data.standard$gene1))
  gene.indices.std <- out 
  unique.names.genes.std <- row.names(out)
  
  # Now sort the similarity data by gene pairs
  if (is.unsorted(data.interaction[,1])){
    ind <- order(data.interaction[,1])
    data.interaction <- data.interaction[ind,]
  }
  out <- GroupUniqueElements(data.interaction[,1])
  gene.indices.sim <- out 
  unique.names.genes.sim <- row.names(out)
  
  
  # Fill out similarity values and corresponding co-annotation (1/0) for each
  # common gene pairs between similarity and co-annotation
  combined.true <- numeric(dim(data.standard)[1]) # 0 - default/ initialization value
  combined.score <- numeric(dim(data.standard)[1])
  source <- vector('character', dim(data.standard)[1])
  indices.in.standard <- numeric(dim(data.standard)[1])
  
  ## Speedup: Instead of using rownames directly, find a numeric representation
  int.genes.in.std <- intersect(unique.names.genes.sim, unique.names.genes.std)
  ind.int.genes.in.std <- match(int.genes.in.std, unique.names.genes.sim)
  
  curr.ind <- 1 # Track the progress
  
  # For pairwise corr, the last one won't be there
  print('Associating similarity values of pairs to pos(1) / neg(0) co-annotation ... ')
  pb <- txtProgressBar(style = 3) # Progress bar
  for (i in ind.int.genes.in.std){
    # Query gene entrez ID
    # common.gene.ID <- gene.symbol[i] # gene symbol
    # curr.unique.gene <- intersect(common.gene.ID, unique.names.genes.std)
    
    curr.unique.gene <- unique.names.genes.sim[i]
    
    ## Get co-annotation values for all pairs of curr.unique.gene
    ind.curr <- gene.indices.std[curr.unique.gene,1] : gene.indices.std[curr.unique.gene,2]
    std.second.genes <- data.standard[ind.curr, 2]
    co.ann.values <- data.standard[ind.curr, 3]
    names(co.ann.values) <- std.second.genes
    
    ## Get similarity data for all pairs of curr.unique.gene
    ind.sim <- gene.indices.sim[i,1] : gene.indices.sim[i,2]
    similarity.values <- data.interaction[ind.sim, 3] # parentheses is important
    names(similarity.values) <- data.interaction[ind.sim, 2]
    
    # Remove the nan correlation values and corresponding interactions
    ind.nan <- which(is.nan(similarity.values)) # Get the indices of NA's
    if (length(ind.nan) > 0){
      similarity.values <- similarity.values[-ind.nan]
    }
    
    ## Assign the true lables (1,0) from the functional standard
    common.genes <- intersect(names(similarity.values), std.second.genes)
    values.true <- co.ann.values[common.genes]
    values.predicted <- similarity.values[common.genes]
    
    curr.size <- length(values.true)
    combined.true[curr.ind: (curr.ind + curr.size - 1)]  <-  values.true
    combined.score[curr.ind: (curr.ind + curr.size - 1)] <-  values.predicted
    
    # If we have to store the IDs and indices too
    if (dim(data.standard)[2] == 4){
      co.ann.IDs <- data.standard[ind.curr, 4]
      names(co.ann.IDs) <- std.second.genes
      source[curr.ind: (curr.ind + curr.size - 1)] <- co.ann.IDs[common.genes]
      
      values.indices <- ind.curr
      names(values.indices) <- std.second.genes
      indices.in.standard[curr.ind: (curr.ind + curr.size - 1)] <- values.indices[common.genes]
    }
    
    curr.ind <- curr.ind + curr.size # Update
    
    setTxtProgressBar(pb, i/length(ind.int.genes.in.std)) # Progress bar update
  }
  close(pb)
  
  ## Remove the unnecessary part of the data
  combined.true <- combined.true[1: (curr.ind - 1)]
  combined.score <- combined.score[1: (curr.ind - 1)]
  
  ## Remove the na correlation values and corresponding interactions
  ind.na <- which(is.na(combined.score) | is.nan(combined.score))
  if (length(ind.na) > 0){
    combined.score <- combined.score[-ind.na]
    combined.true <- combined.true[-ind.na]  
  }
  
  # If data.standard is of dimension 4 (contains ID information)
  if (dim(data.standard)[2] == 4){
    source <- source[1: (curr.ind - 1)]
    indices.in.standard <- indices.in.standard[1: (curr.ind - 1)]
    
    if (length(ind.na) > 0){
      source <- source[-ind.na]
      indices.in.standard <- indices.in.standard[-ind.na]
    }
  }
  
  # Returning as a list
  if (dim(data.standard)[2] == 4){
    return(list(true = combined.true, predicted = combined.score, ID = source, index = indices.in.standard))  
  } else{
    return(list(true = combined.true, predicted = combined.score))
  }
}



#' Pair up the direct interaction between different pairs and their true co-annotation value
#'
#' @param data.standard: the functional standard (Symbol1, Symbol2, 1/0, source (optional))
#' @param data.interaction: interaction scores (returned from GetInteractionData)
#'                   or pairwise.correlation on interaction data
#' @return 
#'  A list of :
#'   data: combined (query) interactions (score) and their corresponding co-annotation (true) (sorted from -ve to +ve)
#'   Queries.with.AUC: The queries that has some AUC values (others are 0), sorted alphabetically
#'   auc.neg: AUC value for the data (neg to pos), sorted according to queries
#'   auc.pos: AUC value for the data (pos to neg), sorted according to queries
#' @export
#' 
CalculatePredictionAndTrueOnDirectInteraction <- function(data.standard, data.interaction){

  ## Subset the data using unique query genes (if not already provided) 
  queries <- unlist(lapply(strsplit(colnames(data.interaction), '_'), '[[', 1) ) # Get the first element of the list (strings separated by '_')
  queries_unique <- unique(queries)

  # If there are multiple replicates, we are taking the first unique one  
  indices <- match(queries_unique, queries)
  data.interaction <- data.interaction[,indices]
  queries <- queries[indices]
  queries.original <- colnames(data.interaction)
  colnames(data.interaction) <- queries # queries.original[indices]
  queries.full <- queries.original[indices] # Complete names of the query

  ## PR Curve for each query screen (using numeric interaction data)
  combined.true.neg <- c()
  combined.score.neg <- c()
  GI_ones <- numeric(length(queries)) # Number of interactions for a gene
  AUC.neg.summary <- c()
  AUC.pos.summary <- c()
  Queries.with.AUC <- c()
  
  
  # Get the indices of unique genes (sorted by the first gene of the pair)
  if (is.unsorted(data.standard$gene1)){
    ind <- order(data.standard$gene1)
    data.standard <- data.standard[ind,]
  }
  
  # system.time(out <- GroupUniqueElements(data.standard$gene1))
  out <- GroupUniqueElements(data.standard$gene1)
  indices.genes1 <- out 
  unique.genes1 <- row.names(out)
  
  # This is for the incidences where the gene name is on gene2
  ind <- order(data.standard$gene2)
  data.standard.tmp <- data.standard[ind,]
  
  out <- GroupUniqueElements(data.standard.tmp$gene2)
  indices.genes2 <- out 
  unique.genes2 <- row.names(out)
  
  curr.ind <- 0
  
  ## For all queries
  pb <- txtProgressBar(style = 3) # Progress bar
  for (query.cand in queries){
    # query.cand = 'MUS81' # Complex testing
    # query.cand = 'ABCC1' # Pathway testing
    
    curr.ind <- curr.ind + 1
    
    # If the query is absent from the standard, skip it
    if (!(is.element(query.cand, unique.genes1) | is.element(query.cand, unique.genes2))){
      next
    }
    else{ # Get the standard data for this query Gene
      # print(query.cand)
      
      func.inter.genes <- c()
      func.inter.values <- c()
      
      # If it comes from first gene (data.standard)
      if (is.element(query.cand, unique.genes1)){
        ind1 <- indices.genes1[query.cand,1] : indices.genes1[query.cand,2]  
        func.inter.genes <- data.standard[ind1,]$gene2 # unique(data.standard[ind1,1])
        func.inter.values <- data.standard[ind1,]$is_annotated
      }
      
      # If it comes from second gene (data.standard.tmp)
      if (is.element(query.cand, unique.genes2)){
        ind2 <- indices.genes2[query.cand,1] : indices.genes2[query.cand,2]  
        func.inter.genes <- c(func.inter.genes, data.standard.tmp[ind2,]$gene1) # unique(data.standard.tmp[ind2,2])
        func.inter.values <- c(func.inter.values, data.standard.tmp[ind2,]$is_annotated)
      }
      
      names(func.inter.values) <- func.inter.genes
      GI_ones <- c(GI_ones, sum(func.inter.values == 1))
    }
    
    genetic.inter.values <- data.interaction[,query.cand] # Single column
    names(genetic.inter.values) <- row.names(data.interaction)
    
    ## Pair interactions / scores with co-annotations
    common.pairs <- intersect(names(genetic.inter.values), names(func.inter.values))
    
    tmp.true.neg <- func.inter.values[common.pairs]
    tmp.score.neg <- genetic.inter.values[common.pairs]
    
    if( length(unique(tmp.true.neg)) == 2) { # we have both classes
      
      # Sort (ascendingly): Neg to Pos (Not required actually)
      ind.ascending <- order(tmp.score.neg)
      tmp.score.neg <- tmp.score.neg[ind.ascending]
      tmp.true.neg <- tmp.true.neg[ind.ascending]
      
      # Save globally for combined query analysis later
      combined.true.neg  = c(combined.true.neg, tmp.true.neg)
      combined.score.neg = c(combined.score.neg, tmp.score.neg)
      
      # AUC (neg to pos): ascend = TRUE
      output.neg <- GenerateDataForPerfCurve(value.predicted = tmp.score.neg, value.true = tmp.true.neg, 
                                             neg.to.pos = TRUE, x.axis = 'recall', y.axis = 'precision')
      AUC.neg.summary <- c(AUC.neg.summary, output.neg$auc)
      
      # AUC (pos to neg): ascend = FALSE
      output.pos <- GenerateDataForPerfCurve(value.predicted = tmp.score.neg, value.true = tmp.true.neg, 
                                             neg.to.pos = FALSE, x.axis = 'recall', y.axis = 'precision')
      AUC.pos.summary <- c(AUC.pos.summary, output.pos$auc)
      
      # Queries.with.AUC <- c(Queries.with.AUC, query.cand)
      Queries.with.AUC <- c(Queries.with.AUC, queries.full[curr.ind])
      print(curr.ind)
    }
    else{ # We just have one class (so no PR curve)
      next 
    }
    
    setTxtProgressBar(pb, curr.ind/length(query.cand)) # Progress bar update
  }
  close(pb)
  
  ## Remove NaN from the data
  ind.na <- which(is.na(combined.score.neg) | is.nan(combined.score.neg))
  if (length(ind.na) > 0){
    combined.score.neg <- combined.score.neg[-ind.na]
    combined.true.neg = combined.true.neg[-ind.na]
  }
  
  ## Order from negative to positive score
  tmpInd = order(combined.score.neg)
  combined.score.neg <- combined.score.neg[tmpInd]
  combined.score.neg = combined.score.neg[tmpInd]
  combined.neg = data.frame(True = combined.true.neg, Score = combined.score.neg)
  
  ## Sort by Query Names
  ind <- order(Queries.with.AUC)
  Queries.with.AUC <- Queries.with.AUC[ind]
  AUC.neg.summary <- AUC.neg.summary[ind]
  AUC.pos.summary <- AUC.pos.summary[ind]
  
  return(list(data = combined.neg, queries = Queries.with.AUC, auc.neg.query = AUC.neg.summary, auc.pos.query = AUC.pos.summary))
}
