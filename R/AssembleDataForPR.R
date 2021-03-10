## Functions to associate pairwise scores with ground truth values.
## Copyright (C) 2018-2021 AHM Mahfuzur Rahman (rahma118@umn.edu)
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
  
  # ------------- 1. data.standard -------------
  true_classes_count <- ((class(data.standard) == 'data.frame') | (class(data.standard) == 'list')) + (class(data.interaction) == 'data.frame' | class(data.interaction) == 'matrix')
  if(true_classes_count < 2){
    stop ("Check the inputs again ... they are supposed to be data.frame's")
  }
  
  true_column_type <- ((class(data.standard[,1]) == 'character') + (class(data.standard[,2]) == 'character') + ((class(data.standard[,3]) == 'integer') | (class(data.standard[,3]) == 'numeric')))
  if(true_column_type < 3){
    stop ("Check the column types. They should be 'character', 'character', 'integer'. Fourth column is optional ...")
  }
  
  # ------------- 2. data.interaction -------------
  # Format: rownames are genes and the actual values are numbers
  
  # Check if all columns are numeric
  if(!sum(unlist(lapply(data.interaction, is.numeric))) == dim(data.interaction)[2]){
    stop ("data.interaction contains non-numeric columns (values)")
  }
  
  if (is.null(row.names(data.interaction[1,,drop=FALSE]))){
    stop ("data.interaction does not contain gene names as row names ...")
  }
  
  ## *** Decide which function to call and prepare data for that!
  min_dim <- min(dim(data.interaction)[1], dim(data.interaction)[2])
  subset <- sort(sample(1:min_dim, min_dim / 10))
  
  # *** 1. Pairwise data in the format of Gene1 Gene2 Similarity
  if (dim(data.interaction)[2] == 3) {
    # We want to check the data types for this format
    if ((typeof(data.interaction[1,1]) == 'character') & 
        (typeof(data.interaction[1,2]) == 'character') & 
        (typeof(data.interaction[1,3]) == 'double')){
      
      print('Similarity data provided ... calling FromGenePairSimilarity ...')
      
      # Use Symbol or Entrez ID (for GIANT)
      if (class(data.standard) == 'list'){
        return (FromGenePairSimilarityEntrez(data.standard, data.interaction))
      } else{
        return (FromGenePairSimilarity(data.standard, data.interaction))
      }
      
    } else{
      print('ill-formatted data .. quitting ...')
      return(NULL)
    }
    
  } else if ( (dim(data.interaction)[1] == dim(data.interaction)[2]) & 
              isSymmetric(as.matrix(data.interaction[subset, subset]))) {
    # *** 2. Pairwise correlation / similarity data provided as a square matrix    
    # Check: 1. Square, 2. Symmetric (We are using a random subset to reduce memory footprint)
    print('Pairwise correlation/similarity matrix provided ...')
    
    # Sort the pairwise data by gene names (on both rows and columns)
    if (is.unsorted(rownames(data.interaction))){
      gene.names <- rownames(data.interaction)
      ind <- order(gene.names)
      data.interaction <- data.interaction[ind,ind]      
    }
    data.interaction <- as.matrix(data.interaction) # This is pairwise.corr
    
  } else{ 
    # *** 3. Genes (row) * screens (column) matrix provided (GI or fitness, for example)
    # Calculate pairwise.correlation
    
    # Sort the interaction.data by gene names (library side)
    if (is.unsorted(rownames(data.interaction))){
      gene.names <- rownames(data.interaction)
      ind <- order(gene.names)
      data.interaction <- data.interaction[ind,]
      row.names(data.interaction) <- gene.names[ind]
    }
    
    print('calculating pairwise correlation ...')
    
    # Check if there are a lot of NaNs in data? Total percentage, or percent of NaN columns
    percent.nan <- sum(is.nan(as.matrix(data.interaction))) / (dim(data.interaction)[1] * dim(data.interaction)[2]) + sum(is.na(as.matrix(data.interaction))) / (dim(data.interaction)[1] * dim(data.interaction)[2])
    percent.nan.cols <- (sum(is.nan(apply(data.interaction, sum, MARGIN = 2))) + sum(is.na(apply(data.interaction, sum, MARGIN = 2)))) / dim(data.interaction)[2]
    
    if ( (percent.nan > 0.1) | (percent.nan.cols > 0.1) ){ # If we have more than 20% NaNs in the data
      # data.interaction <- cor(t(data.interaction), use = 'pairwise.complete.obs', method = 'pearson') # Slow, and really unstable (for sparse data)
      
      ## Impute nan/NA values
      # Way 1: using a random value based on background (May be more variable)
      # int_mean <- mean(as.vector(as.matrix(data.interaction)), na.rm = T) # mean is very close to 0 anyway!
      # int_sd <- sd(as.vector(as.matrix(data.interaction)), na.rm = T)
      # num_na <- sum(is.na(as.matrix(data.interaction)))
      # if(num_na > 0)
      #  data.interaction[is.na(as.matrix(data.interaction))] <- rnorm(num_na, mean = 0, sd = int_sd)
      
      # Way 2: Use a mean/median per gene score?
      k <- which(is.na(data.interaction), arr.ind=TRUE)
      data.interaction[k] <- rowMeans(data.interaction, na.rm=TRUE)[k[,1]]
      
      # Now calcualte correlation: no NA in the set
      data.interaction <- cor(t(data.interaction), use = 'complete.obs', method = 'pearson')
      
    }else{
      data.interaction <- cor(t(data.interaction), use = 'complete.obs', method = 'pearson') # Fast
    }
    # pairwise.correlation <- data.interaction
  }
  
  # Use Symbol or Entrez ID (for GIANT)
  if (class(data.standard) == 'list'){
    return (FromAllPairwiseCorrelationEntrez(data.standard, data.interaction))
  } else{
    return (FromAllPairwiseCorrelation(data.standard, data.interaction))
  }
  
}


# Supporting function for CalculatePredictionAndTrueOnLibraryProfiles
FromAllPairwiseCorrelation <- function(data.standard, pairwise.correlation){
  
  print('In FromAllPairwiseCorrelation (data.standard, data.interaction) ...')
  
  ## Pre-processing ===================================
  # print('speed optimizations ...')
  length.corr <- dim(pairwise.correlation)[1]
  gene.symbol <- rownames(pairwise.correlation)
  
  # Sort the co annotation data (data.standard) if not already sorted
  if (is.unsorted(data.standard$gene1)){
    ind <- order(data.standard$gene1)
    data.standard <- data.standard[ind,]
  }
  
  gene.indices <- GroupUniqueElements(data.standard$gene1)
  unique.names.genes <- row.names(gene.indices)
  
  # Fill out similarity values and corresponding co-annotation (1/0) for each
  # common gene pairs between similarity and co-annotation
  combined.true <- numeric(dim(data.standard)[1]) # 0 - default/ initialization value
  combined.score <- numeric(dim(data.standard)[1])
  source <- vector('character', dim(data.standard)[1])
  indices.in.standard <- numeric(dim(data.standard)[1])
  
  ## Probable Speedup: Not use rownames directly, but find some numeric representation
  int.genes.in.std <- intersect(gene.symbol, unique.names.genes)
  ind.int.genes.in.std <- match(int.genes.in.std, gene.symbol) # 17144
  # identical(gene.symbol[ind.int.genes.in.std], int.genes.in.std) # TRUE
  
  
  ## Main Processing ==============================================
  #  Loop through all library genes that are also in the standard
  count <- 0
  
  curr.ind <- 1 # Track the progress
  print('Mapping pairs against their correlation values ... ')
  
  # Progress bar
  pb <- txtProgressBar(style = 3) 
  pb_count <- 0
  
  for (i in ind.int.genes.in.std){ # i -> index of gene in data.interaction
    curr.unique.gene <- gene.symbol[i]
    
    ## *** The last row doesn't matter, as we have already covered the gene.
    if (i == length.corr){
      next # This basically means break in this case!
    }
    
    ## Get standard data for this names.genes (with others associated with it)
    ind.pairs.in.std <- gene.indices[curr.unique.gene,1] : gene.indices[curr.unique.gene,2]
    
    std.second.genes <- data.standard[ind.pairs.in.std, 2] # gene 2
    co.ann.values <- data.standard[ind.pairs.in.std, 3] # co-ann
    names(co.ann.values) <- std.second.genes
    values.indices <- ind.pairs.in.std # index as a pointer to original standard
    names(values.indices) <- std.second.genes
    
    # *** Way1 : Considering only Upper triangle (both datasets must be sorted)
    # data.corr <- pairwise.correlation[i, (i+1) : length.corr] # parentheses is important
    # sim.second.genes <- gene.symbol[(i+1) : length.corr]
    # print(length(intersect(sim.second.genes, std.second.genes)))
    
    # *** Way2: Consdering whole
    # Because for standards not created by us, like GIANT, we may lose some pairs! 
    # Problem: If the standard has duplicate pairs, they will affect the PR performance.
    data.corr <- pairwise.correlation[i, ]
    sim.second.genes <- names(data.corr) # This works too!
    # print(length(intersect(sim.second.genes, std.second.genes)))
    
    # system.time(sim.second.genes <- gene.symbol[(i+1) : length.corr]) # Faster
    # system.time(sim.second.genes <- gene.symbol[-(1:i)])
    
    ## Assign the true lables (1,0) from the functional standard
    common.genes <- intersect(sim.second.genes, std.second.genes)
    values.true <- co.ann.values[common.genes]
    values.predicted <- data.corr[common.genes]
    values.indices <- values.indices[common.genes]
    curr.size <- length(values.true)
    
    # print(curr.size)
    # print(paste(i, length(values.true), sep = ' '))
    
    if (curr.size > 0){ # R will produce an error otherwise
      combined.true[curr.ind : (curr.ind + curr.size - 1)]  <-  values.true
      combined.score[curr.ind : (curr.ind + curr.size - 1)] <-  values.predicted
      indices.in.standard[curr.ind: (curr.ind + curr.size - 1)] <- values.indices
      
      # If the ID of the pair (which complex/PW it belongs) is provided
      if (dim(data.standard)[2] == 4){
        co.ann.IDs <- data.standard[ind.pairs.in.std, 4]
        names(co.ann.IDs) <- std.second.genes
        source[curr.ind: (curr.ind + curr.size - 1)] <- co.ann.IDs[common.genes]
      }
      count <- count + 1
    }
    
    curr.ind <- curr.ind + curr.size # Update
    
    # Progress bar update
    pb_count <- pb_count + 1
    setTxtProgressBar(pb, pb_count/length(ind.int.genes.in.std))
  }
  close(pb)
  
  ## *** Logical error check
  if (curr.ind < 2){
    return (NULL)
  }
  
  ## Post-processing ===================================
  ## Remove the unnecessary part of the data
  combined.true <- combined.true[1: (curr.ind - 1)] # 7220816 (GIANT) (the other way gives 7235130)
  combined.score <- combined.score[1: (curr.ind - 1)]
  indices.in.standard <- indices.in.standard[1: (curr.ind - 1)]
  
  ## Remove the na correlation values and corresponding interactions
  ind.na <- which(is.na(combined.score) | is.nan(combined.score))
  if (length(ind.na) > 0){
    combined.score <- combined.score[-ind.na]
    combined.true <- combined.true[-ind.na]
    indices.in.standard <- indices.in.standard[-ind.na]
  }
  
  # If data.standard is of dimension 4 (contains ID information)
  if (dim(data.standard)[2] == 4){
    source <- source[1: (curr.ind - 1)]
    if (length(ind.na) > 0){
      source <- source[-ind.na]
    }
  }
  
  # Sanity check (should be TRUE)
  # identical(combined.true, data.standard[indices.in.standard,3])
  
  # Return outputs as a list
  if (dim(data.standard)[2] == 4){
    return(list(true = combined.true, predicted = combined.score, ID = source, index = indices.in.standard))
  } else{
    return(list(true = combined.true, predicted = combined.score, index = indices.in.standard))
  }
  
}


# Supporting function for CalculatePredictionAndTrueOnLibraryProfiles
FromGenePairSimilarity <- function(data.standard, data.interaction){
  
  print('In FromGenePairSimilarity (data.standard, data.interaction) ...')
  
  ## Pre-processing ===================================
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
  
  ## Main Processing ==============================================
  #  Loop through all library genes that are also in the standard
  for (i in ind.int.genes.in.std){
    curr.unique.gene <- unique.names.genes.sim[i]
    
    ## *** The last row doesn't matter, as we have already covered the gene.
    # if (i == length.corr){
    #  next # This basically means break in this case!
    #}
    
    ## Get co-annotation values for all pairs of curr.unique.gene
    ind.pairs.in.std <- gene.indices.std[curr.unique.gene,1] : gene.indices.std[curr.unique.gene,2]
    
    std.second.genes <- data.standard[ind.pairs.in.std, 2] # gene 2
    co.ann.values <- data.standard[ind.pairs.in.std, 3] # co-ann
    names(co.ann.values) <- std.second.genes
    values.indices <- ind.pairs.in.std # index as a pointer to original standard
    names(values.indices) <- std.second.genes
    
    ## Get similarity data for all pairs of curr.unique.gene
    ind.sim <- gene.indices.sim[i,1] : gene.indices.sim[i,2]
    similarity.values <- data.interaction[ind.sim, 3]
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
    values.indices <- values.indices[common.genes]
    
    curr.size <- length(values.true)
    
    if (curr.size > 0){ # R will produce an error otherwise
      
      combined.true[curr.ind: (curr.ind + curr.size - 1)]  <-  values.true
      combined.score[curr.ind: (curr.ind + curr.size - 1)] <-  values.predicted
      indices.in.standard[curr.ind: (curr.ind + curr.size - 1)] <- values.indices
      
      # If the ID of the pair (which complex/PW it belongs) is provided
      if (dim(data.standard)[2] == 4){
        co.ann.IDs <- data.standard[ind.pairs.in.std, 4]
        names(co.ann.IDs) <- std.second.genes
        source[curr.ind: (curr.ind + curr.size - 1)] <- co.ann.IDs[common.genes]
      }
    }    
    
    curr.ind <- curr.ind + curr.size # Update
    setTxtProgressBar(pb, i/length(ind.int.genes.in.std)) # Progress bar update
  }
  close(pb)
  
  ## *** Logical error check
  if (curr.ind < 2){
    return (NULL)
  }
  
  ## Post-processing ===================================
  ## Remove the unnecessary part of the data
  combined.true <- combined.true[1: (curr.ind - 1)]
  combined.score <- combined.score[1: (curr.ind - 1)]
  indices.in.standard <- indices.in.standard[1: (curr.ind - 1)]
  
  ## Remove the na correlation values and corresponding interactions
  ind.na <- which(is.na(combined.score) | is.nan(combined.score))
  if (length(ind.na) > 0){
    combined.score <- combined.score[-ind.na]
    combined.true <- combined.true[-ind.na]
    indices.in.standard <- indices.in.standard[-ind.na]
  }
  
  # If data.standard is of dimension 4 (contains ID information)
  if (dim(data.standard)[2] == 4){
    source <- source[1: (curr.ind - 1)]
    if (length(ind.na) > 0){
      source <- source[-ind.na]
    }
  }
  
  # Returning as a list
  if (dim(data.standard)[2] == 4){
    return(list(true = combined.true, predicted = combined.score, ID = source, index = indices.in.standard))  
  } else{
    return(list(true = combined.true, predicted = combined.score, index = indices.in.standard))
  }
}



#' Pair up the direct interaction between different gene pairs and their true co-annotation value
#'
#' @param data.standard: the functional standard (Symbol1, Symbol2, 1/0, source (optional))
#' @param data.interaction: interaction scores (returned from GetInteractionData)
#' @param provide.identity: if TRUE, return gene pairs for each interaction (library and query)
#'
#' @return 
#'  A list of :
#'   data: combined (query) interactions (score) and their corresponding co-annotation (true) (sorted from -ve to +ve)
#'   Queries.with.AUC: The queries that has some AUC values (others are 0), sorted alphabetically
#'   auc.neg: AUC value for the data (neg to pos), sorted according to queries
#'   auc.pos: AUC value for the data (pos to neg), sorted according to queries
#' @export
#' 
CalculatePredictionAndTrueOnDirectInteraction <- function(data.standard, data.interaction, provide.identity = FALSE){
  
  ## 1. Sort and group data from standard
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
  
  
  ## 2. Subset the data using unique query genes (if not already provided) 
  queries <- unlist(lapply(strsplit(colnames(data.interaction), '_'), '[[', 1) ) # Get the first element of the list (strings separated by '_')
  queries_unique <- unique(queries)
  
  # If there are multiple replicates, we are taking the first unique one  
  indices <- match(queries_unique, queries)
  data.interaction <- data.interaction[,indices]
  queries <- queries[indices]
  queries.original <- colnames(data.interaction)
  colnames(data.interaction) <- queries # queries.original[indices]
  queries.full <- queries.original[indices] # Complete names of the query
  
  ## 3. PR Curve for each query screen (using numeric interaction data)
  combined.true.neg <- c()
  combined.score.neg <- c()
  
  if(provide.identity == TRUE){
    combined.library.gene <- c()
    combined.query.gene <- c()
  }
  
  GI_ones <- c() # Number of interactions for a query
  AUC.neg.summary <- c()
  AUC.pos.summary <- c()
  Queries.with.AUC <- c()
  
  ## 4. Generate the Score vs Annotation data (all queries)
  curr.ind <- 0
  pb <- txtProgressBar(style = 3) # Progress bar
  for (query.cand in queries){
    # query.cand = 'MUS81' # Complex testing
    # query.cand = 'ABCC1' # Pathway testing
    # query.cand = 'AKT1'
    
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
      combined.true.neg  <- c(combined.true.neg, tmp.true.neg)
      combined.score.neg <- c(combined.score.neg, tmp.score.neg)
      
      if(provide.identity == TRUE){
        # Store gene names
        combined.library.gene <- c(combined.library.gene, names(tmp.true.neg))
        combined.query.gene <- c(combined.query.gene, rep(query.cand, length(tmp.true.neg)))
      }
      
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
    combined.true.neg <- combined.true.neg[-ind.na]
  }
  
  ## Order from negative to positive score
  tmpInd = order(combined.score.neg)
  combined.score.neg <- combined.score.neg[tmpInd]
  combined.true.neg <- combined.true.neg[tmpInd]
  
  if(provide.identity == TRUE){
    # Store gene names
    combined.library.gene <- combined.library.gene[tmpInd]
    combined.query.gene <- combined.query.gene[tmpInd]
    combined.neg = data.frame(True = combined.true.neg, Score = combined.score.neg, 
                              library = combined.library.gene, query = combined.query.gene)
  } else{
    # Don't store gene names
    combined.neg = data.frame(True = combined.true.neg, Score = combined.score.neg)
  }
  
  ## Sort by Query Names
  ind <- order(Queries.with.AUC)
  Queries.with.AUC <- Queries.with.AUC[ind]
  AUC.neg.summary <- AUC.neg.summary[ind]
  AUC.pos.summary <- AUC.pos.summary[ind]
  
  return(list(data = combined.neg, queries = Queries.with.AUC, auc.neg.query = AUC.neg.summary, auc.pos.query = AUC.pos.summary))
}


## ==============================================================================
##          Standard in Entrez, interaction in Symbol
## ==============================================================================
# Supporting function for CalculatePredictionAndTrueOnLibraryProfiles
FromAllPairwiseCorrelationEntrez <- function(data.standard, pairwise.correlation){
  
  print('In FromAllPairwiseCorrelationEntrez (data.standard, data.interaction) ...')
  
  ## Pre-processing ===================================
  if (class(data.standard) == 'list'){
    ## Get the groupings
    gene.indices <- data.standard$gene.indices # sorted by entrez ID
    unique.genes.std.entrez <- row.names(gene.indices)
    
    ## Get the mappings
    mapping <- data.standard$mapping
    unique.genes.std.symbol <- mapping[row.names(gene.indices)]
    
    ## Main data
    data.ca <- data.standard$data
  } else{
    stop('Expects mapping information ....')
  }
  
  # Convert the Gene symbols to Entrez
  # length.corr <- dim(pairwise.correlation)[1]
  common.genes <- intersect(rownames(pairwise.correlation), mapping) 
  pairwise.correlation <- pairwise.correlation[common.genes, common.genes]
  gene.symbol <- rownames(pairwise.correlation) # 17206 (out of 17634)
  
  ## Speedup: Not use rownames directly, but find some numeric representation
  common.genes <- intersect(gene.symbol, unique.genes.std.symbol) # 17148
  ind.int.genes <- match(common.genes, gene.symbol)
  ind.std.genes <- match(common.genes, unique.genes.std.symbol)
  # points to genes alphabetically "A1BG" "A1CF" "A2M" "A2ML1" etc.
  
  ## Main Processing ==============================================
  print('Mapping pairs against their correlation values ... ')
  count <- 0
  
  combined.true <- numeric(dim(data.ca)[1]) 
  combined.score <- numeric(dim(data.ca)[1])
  indices.in.standard <- numeric(dim(data.ca)[1])
  
  curr.ind <- 1 # Track the progress
  
  pb <- txtProgressBar(style = 3) # Progress bar
  pb_count <- 0
  
  for (i in seq_along(ind.int.genes)){ # i -> index of gene in data.interaction
    
    # GREAT !!! Getting values for same genes from both
    # i <- seq_along(ind.int.genes)
    # ind.int.genes[i] # gene index here
    # gene.symbol[ind.int.genes[i]] # gene symbols
    # ind.std.genes[i] # entrez ID
    # mapping[row.names(gene.indices)[ind.std.genes[i]]]
    # identical(gene.symbol[ind.int.genes[i]], unname(mapping[row.names(gene.indices)[ind.std.genes[i]]])) # TRUE
    
    # i <- 3 # 'A2M'
    
    curr.int.gene <- ind.int.genes[i] # apply as indices to pairwise.correlation
    curr.std.gene <- ind.std.genes[i] # appy as indices to gene.indices
    
    ## Get genes associated to curr.std.gene
    ind.pairs.in.std <- gene.indices[curr.std.gene,1] : gene.indices[curr.std.gene,2] # unique(data.ca[ind.pairs.in.std,1])
    std.second.genes <- mapping[as.character(data.ca[ind.pairs.in.std, 2])]
    
    non.na.ind <- !is.na(std.second.genes) # There maybe some NA's (no mapping)
    std.second.genes <- std.second.genes[non.na.ind]
    # sum(non.na.ind) # length(non.na.ind)
    
    values.indices <- ind.pairs.in.std[non.na.ind] # Indices of non-NA data from std
    names(values.indices) <- std.second.genes
    # identical(mapping[as.character(data.ca[values.indices,2])], std.second.genes)
    
    co.ann.values <- data.ca[values.indices, 3]
    names(co.ann.values) <- std.second.genes
    
    ## Get the interaction data
    data.corr <- pairwise.correlation[curr.int.gene, ]
    int.second.genes <- names(data.corr)
    
    ## Assign the true lables (1,0) from the functional standard
    common.genes <- intersect(int.second.genes, std.second.genes)
    curr.size <- length(common.genes)
    
    # print(paste(unique.genes.std.entrez[curr.std.gene], curr.size, sep = ' '))
    
    if (curr.size <= 0){ next }
    
    combined.true[curr.ind:(curr.ind+curr.size-1)]  <-  co.ann.values[common.genes] # 1/0
    combined.score[curr.ind:(curr.ind+curr.size-1)] <-  data.corr[common.genes] # Corr
    indices.in.standard[curr.ind:(curr.ind+curr.size-1)] <- values.indices[common.genes]
      
    count <- count + 1
    curr.ind <- curr.ind + curr.size # Update
    
    pb_count <- pb_count + 1 # Progress bar update
    setTxtProgressBar(pb, pb_count/length(ind.int.genes))
  }
  close(pb)
  
  ## *** Logical error check
  if (curr.ind < 2){
    return (NULL)
  }
  
  ## Post-processing ===================================
  combined.true <- combined.true[1: (curr.ind - 1)]
  combined.score <- combined.score[1: (curr.ind - 1)]
  indices.in.standard <- indices.in.standard[1: (curr.ind - 1)]
  
  ## Remove the na correlation values and corresponding interactions
  ind.na <- which(is.na(combined.score) | is.nan(combined.score))
  if (length(ind.na) > 0){
    combined.score <- combined.score[-ind.na]
    combined.true <- combined.true[-ind.na]
    indices.in.standard <- indices.in.standard[-ind.na]
  }
  
  ## Do we have duplicate pairs? (if not indices.in.standard will be unique)
  sum(duplicated(indices.in.standard)) # FALSE :) ... Check again!
  
  # Return outputs as a list
    return(list(true = combined.true, predicted = combined.score, index = indices.in.standard))
  
  # Sanity Check! It's working !!!
  # data.out <- data.frame(true = combined.true, predicted = combined.score, index = indices.in.standard)
  # identical(data.out$true, data.ca[data.out$index, 3]) # TRUE
  # pairwise.correlation[mapping[as.character(data.ca[head(data.out)$index,]$gene1)],mapping[as.character(data.ca[head(data.out)$index,]$gene2)]]
  # pairwise.correlation[mapping[as.character(data.ca[tail(data.out)$index,]$gene1)],mapping[as.character(data.ca[tail(data.out)$index,]$gene2)]]
  
  # Profile_DepMap_19Q2 = list(true = combined.true, predicted = combined.score, index = indices.in.standard)
  # pred.ca <- list(test = list(true = Profile_DepMap_19Q2$true, predicted = Profile_DepMap_19Q2$predicted))
  # PlotPRSimilarity (pred.ca, subsample = TRUE, type.plot = 'log', fig.title = paste0('Profile Similarity (GIANT)'), legend.names = c('19Q2'), legend.color = c('#756bb1'), save.figure = TRUE, outfile.name = paste0('DepMap_19Q2_Profile'),  outfile.type = 'png')
}


# Supporting function for CalculatePredictionAndTrueOnLibraryProfiles
FromGenePairSimilarityEntrez <- function(data.standard, data.interaction){
  
  print('In FromGenePairSimilarityEntrez (data.standard, data.interaction) ...')
  
  ## Pre-processing ===================================
  if (class(data.standard) == 'list'){
    ## Get the groupings (by gene 1)
    gene.indices.std <- data.standard$gene.indices # sorted by entrez ID
    unique.genes.std.entrez <- row.names(gene.indices.std)
    
    ## Get the mappings
    mapping <- data.standard$mapping
    unique.genes.std.symbol <- mapping[unique.genes.std.entrez]
    
    ## Main data
    data.ca <- data.standard$data
  } else{
    stop('Expects mapping information ....')
  }
  
  ## Fill out similarity values and corresponding co-annotation (1/0) for each
  # common gene pairs between similarity and co-annotation
  combined.true <- numeric(2 * dim(data.ca)[1]) # 0 - default/ initialization value
  combined.score <- numeric(2 * dim(data.ca)[1])
  source <- vector('character', 2 * dim(data.ca)[1])
  indices.in.standard <- numeric(2 * dim(data.ca)[1])
  
  curr.ind <- 1 # Track the progress
  
  ## *** Sort the similarity data by gene pairs (group by gene1)
  if (is.unsorted(data.interaction[,1])){
    ind <- order(data.interaction[,1])
    data.interaction <- data.interaction[ind,]
  }
  gene.indices.sim <- GroupUniqueElements(data.interaction[,1])
  unique.names.genes.sim <- row.names(gene.indices.sim)
  
  ## Speedup: Instead of using rownames directly, find a numeric representation
  common.genes <- intersect(unique.names.genes.sim, unique.genes.std.symbol) # 17147
  ind.std.genes <- match(common.genes, unique.genes.std.symbol)
  ind.int.genes <- match(common.genes, unique.names.genes.sim)
  
  # For pairwise corr, the last one won't be there
  print('Associating similarity values of pairs to pos(1) / neg(0) co-annotation ... ')
  pb <- txtProgressBar(style = 3) # Progress bar
  
  ## Main Processing ============================================== 
  #  Loop through all library genes that are also in the standard
  for (i in seq_along(ind.int.genes)){
    
    curr.int.gene <- ind.int.genes[i] # apply as indices to pairwise.correlation
    curr.std.gene <- ind.std.genes[i] # appy as indices to gene.indices.std
    
    ## Get co-annotation values for all pairs of curr.std.gene
    ind.pairs.in.std <- gene.indices.std[curr.std.gene,1] : gene.indices.std[curr.std.gene,2]
    std.second.genes <- mapping[as.character(data.ca[ind.pairs.in.std, 2])]
    
    non.na.ind <- !is.na(std.second.genes) # There maybe some NA's (no mapping)
    std.second.genes <- std.second.genes[non.na.ind]
    # sum(non.na.ind) # length(non.na.ind)
    
    values.indices <- ind.pairs.in.std[non.na.ind] # Indices of non-NA data from std
    names(values.indices) <- std.second.genes
    # identical(mapping[as.character(data.ca[values.indices,2])], std.second.genes)
    
    co.ann.values <- data.ca[values.indices, 3]
    names(co.ann.values) <- std.second.genes
    
    ## Get similarity data for all pairs of curr.int.gene
    ind.sim <- gene.indices.sim[curr.int.gene,1] : gene.indices.sim[curr.int.gene,2]
    similarity.values <- data.interaction[ind.sim, 3]
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
    values.indices <- values.indices[common.genes]
    curr.size <- length(values.true)
    
    if (curr.size > 0){ # R will produce an error otherwise
      combined.true[curr.ind: (curr.ind + curr.size - 1)]  <-  values.true
      combined.score[curr.ind: (curr.ind + curr.size - 1)] <-  values.predicted
      indices.in.standard[curr.ind: (curr.ind + curr.size - 1)] <- values.indices
    }    
    
    curr.ind <- curr.ind + curr.size # Update
    setTxtProgressBar(pb, i/length(ind.int.genes)) # Progress bar update
  }
  close(pb)
  
  sum(duplicated(indices.in.standard[1: (curr.ind - 1)])) # No duplicates
  
  ## Second Part: Sort the standard data by gene pairs (group by gene2)
  #  Now for GIANT, and for data not as pairwise.correlation, there is the chance that we are going to miss some pairs, just because both of them have unique pairs and are sorted in a different way (both on gene1, but different IDs). What we are going to do is group on gene2 (effectively flips gene1 and gene2, and sort by that) to get those extra pairs that we missed first time!
  
  ## Get the groupings (by gene 2) - At this point we know that the standard is a list!
  
  # Remove the pairs that are already accounted for!
  # data.ca.unaccounted <- data.ca[setdiff( (1 : dim(data.ca)[1]), indices.in.standard[1: (curr.ind - 1)]), ]
  
  ind.gene2 <- order(data.ca$gene2)
  gene.indices.std <- GroupUniqueElements(data.ca$gene2[ind.gene2]) # 25657 * 2
  unique.genes.std.entrez <- row.names(gene.indices.std)
  unique.genes.std.symbol <- mapping[unique.genes.std.entrez]
  #  data.ca <- data.standard$data ## Not changing the data itself
  
  ## Speedup: Instead of using rownames directly, find a numeric representation
  common.genes <- intersect(unique.names.genes.sim, unique.genes.std.symbol) # 17147
  ind.std.genes <- match(common.genes, unique.genes.std.symbol)
  ind.int.genes <- match(common.genes, unique.names.genes.sim)
  
  # For pairwise corr, the last one won't be there
  print('Part 2: Group standard by gene2 in the pair ...')
  pb <- txtProgressBar(style = 3) # Progress bar
  
  for (i in seq_along(ind.int.genes)){
    
    curr.int.gene <- ind.int.genes[i] # apply as indices to pairwise.correlation
    curr.std.gene <- ind.std.genes[i] # appy as indices to gene.indices.std
    
    ## Get co-annotation values for all pairs of curr.std.gene
    ind.pairs.in.std <- ind.gene2[gene.indices.std[curr.std.gene,1] : gene.indices.std[curr.std.gene,2]] # As we didn't change the data during the grouping
    std.first.genes <- mapping[as.character(data.ca[ind.pairs.in.std, 1])]
    
    non.na.ind <- !is.na(std.first.genes) # There maybe some NA's (no mapping)
    std.first.genes <- std.first.genes[non.na.ind]
    # sum(non.na.ind) # length(non.na.ind)
    
    values.indices <- ind.pairs.in.std[non.na.ind] # Indices of non-NA data from std
    names(values.indices) <- std.first.genes
    # identical(mapping[as.character(data.ca[values.indices,2])], std.first.genes)
    
    co.ann.values <- data.ca[values.indices, 3]
    names(co.ann.values) <- std.first.genes
    
    ## Get similarity data (on interaction) for all pairs of curr.int.gene
    ind.sim <- gene.indices.sim[curr.int.gene,1] : gene.indices.sim[curr.int.gene,2]
    similarity.values <- data.interaction[ind.sim, 3]
    names(similarity.values) <- data.interaction[ind.sim, 2]
    
    # Remove the nan correlation values and corresponding interactions
    ind.nan <- which(is.nan(similarity.values)) # Get the indices of NA's
    if (length(ind.nan) > 0){
      similarity.values <- similarity.values[-ind.nan]
    }
    
    ## Assign the true lables (1,0) from the functional standard
    common.genes <- intersect(names(similarity.values), std.first.genes)
    values.true <- co.ann.values[common.genes]
    values.predicted <- similarity.values[common.genes]
    values.indices <- values.indices[common.genes]
    curr.size <- length(values.true)
    
    if (curr.size > 0){ # R will produce an error otherwise
      combined.true[curr.ind: (curr.ind + curr.size - 1)]  <-  values.true
      combined.score[curr.ind: (curr.ind + curr.size - 1)] <-  values.predicted
      indices.in.standard[curr.ind: (curr.ind + curr.size - 1)] <- values.indices
    }    
    
    curr.ind <- curr.ind + curr.size # Update
    setTxtProgressBar(pb, i/length(ind.int.genes)) # Progress bar update
  }
  close(pb)
  
  ## *** Logical error check
  if (curr.ind < 2){
    return (NULL)
  }  
  
  ## Post-processing ===================================
  ## Remove the unnecessary part of the data
  combined.true <- combined.true[1: (curr.ind - 1)]
  combined.score <- combined.score[1: (curr.ind - 1)]
  indices.in.standard <- indices.in.standard[1: (curr.ind - 1)]
  
  ## Remove the na correlation values and corresponding interactions
  ind.na <- which(is.na(combined.score) | is.nan(combined.score))
  if (length(ind.na) > 0){
    combined.score <- combined.score[-ind.na]
    combined.true <- combined.true[-ind.na]
    indices.in.standard <- indices.in.standard[-ind.na]
  }
  
  ## Do we have duplicate pairs? Nope! Good news!
  sum(duplicated(indices.in.standard))
  
  # combined.true <- combined.true[!duplicated(indices.in.standard)]
  # combined.score <- combined.score[!duplicated(indices.in.standard)]
  # indices.in.standard <- indices.in.standard[!duplicated(indices.in.standard)]
  
  # Returning as a list
  return(list(true = combined.true, predicted = combined.score, index = indices.in.standard))
  
}
