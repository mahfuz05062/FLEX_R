## Data generation for additional visualization of the output
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


#' Separate the contribution of signal into individual entries (a complex, pathway, GO Biological Process etc.)
#' Works on the result of CalculatePredictionAndTrueOnLibraryProfiles
#'
#' @param entity.matrix -> True, Predicted values with pair association to an entity
#' @param top.numbers -> How much of the signal to consider (from the top)
#' @param summary.standard -> Identity of each entity ('ID'), 'Name', Genes insdie ('Gene'), 'Length'
#' @param data.standard -> Co annotation standard
#'
#' @return The entities and their contribution to the signal, also the top top.number of TPs
#' 
GetContributionOfEntitiesAtCertainTP <- function (entity.matrix, top.number, summary.standard, data.standard){

  # Sort by predicted scores
  sort.ind = order(entity.matrix[,'predicted'], decreasing = TRUE)
  entity.matrix = entity.matrix[sort.ind,]

  # Now find the top TP Pairs from the data
  count.tp.pairs = cumsum(entity.matrix[,'true'])
  min.ind = min(which (count.tp.pairs == top.number))

  tmp.data = entity.matrix[1:min.ind, ] # Get the data until we hit top.number of TP
  tmp.data = tmp.data[tmp.data[,'true'] == 1, ] # Keep only the TP Pairs (with true = 1)

  # Split the source by ';' and get the IDs in lists
  tmp.names = strsplit(tmp.data[,'ID'], ';')
  all.IDs = as.integer(unlist(tmp.names))
  
  # unique.IDs = unique(all.IDs)
  # count.IDs = tabulate(match(all.IDs,  unique.IDs)) # May need to change if data size is big (Not efficient)
  
  # Sort and group the IDs (to find how many times an ID occurs)
  all.IDs <- sort(all.IDs)
  indices.unique.IDs <- GroupUniqueElements(all.IDs)
  
  # Get the complexes and number of pair contribution
  complex.IDs <- as.integer(row.names(indices.unique.IDs))
  count.IDs <- indices.unique.IDs[,2] - indices.unique.IDs[,1] + 1
  row.names(summary.standard) <- summary.standard$ID
  
  contribution.entities = data.frame(ID = complex.IDs, top_pairs = count.IDs, Name = summary.standard[row.names(indices.unique.IDs),2])

  # Sort by contribution count
  ind <- order(contribution.entities[,"top_pairs"], decreasing = TRUE)
  contribution.entities <- contribution.entities[ind,]
  
  # Get the top pairs themselves
  top_pairs <- data.standard[tmp.data$index,]
  
  return(list(contribution.entities, top_pairs))
}


#' Separate the contribution of signal into individual entries (a complex, pathway, GO Biological Process etc.)
#'
#' @param summary.standard -> Identity of each entity ('ID'), 'Name', Genes insdie ('Gene'), 'Length'
#' @param data.standard -> Original Co-annotation data for each entity
#' @param entity.matrix -> Output of CalculatePredictionAndTrueOnLibraryProfiles (true, predicted, source, indices)
#'
#' @return Area under the curve for each unique entity (protein complex for example)
#' 
#' @export
#'   
GetAreaUnderPRCurveForEntities <- function (summary.standard, data.standard, entity.matrix){
  
  ## *** Check input data format
  true_classes_count <- (class(summary.standard) == 'data.frame') + (class(data.standard) == 'data.frame') + (class(entity.matrix) == 'data.frame')
  if(true_classes_count < 3){
    stop ("Check the inputs again ... they are supposed to be data.frame's")
    # stopifnot(true_classes_count < 3)
    # print ("Check the inputs again ... they are supposed to be data.frame's")
    # return (NULL)
  }
  
  # Find the unique Entities that are represented
  tmp.ID <- entity.matrix$ID[entity.matrix$true == 1]
  tmp.names = strsplit(tmp.ID, ';')
  unique.complex.ID = sort(unique(as.integer(unlist(tmp.names))))
  
  if (is.unsorted(data.standard$gene1)){
    ind <- order(data.standard$gene1)
    data.standard <- data.standard[ind,]
  }
  
  data.input = data.standard[entity.matrix$index, ] # Part of the standard that is relevant here (we want to matchup everything with this)
  
  # Group by gene1
  ind_1 <- sort(data.input$gene1, index.return = TRUE)
  data_tmp_1 <- data.input[ind_1$ix,];
  indices_1 <- GroupUniqueElements(data_tmp_1$gene1)
  genes_1 <- row.names(indices_1)
  
  # Group by gene2
  ind_2 <- sort(data.input$gene2, index.return = TRUE)
  data_tmp_2 <- data.input[ind_2$ix,];
  indices_2 <- GroupUniqueElements(data_tmp_2$gene2)
  genes_2 <- row.names(indices_2)
  
  # Pre-processing
  AUC.values = numeric(length(unique.complex.ID))
  common.complex.ID <- intersect(unique.complex.ID, summary.standard$ID)
  matched.ind <- match(common.complex.ID, summary.standard$ID, nomatch = 0) # identical(summary.standard$ID[matched.ind], common.complex.ID)
  
  diff.AUC <- numeric(length(unique.complex.ID)) # background precision / random expectation subtracted
  
  curr.ind <- 1
  
  # require(stringi) # don't need it anymore it seems?!
  # https://rdrr.io/rforge/stringi/
  # sum(stri_detect(data_subset$source, fixed = to.match))
  # data_subset$source[stri_detect(data_subset$source, fixed = to.match)]
  
  pb <- txtProgressBar(style = 3) # Progress bar
  for (i in matched.ind){
    # print(i)
    gene_list = unlist(strsplit(toupper(summary.standard$Genes[i]), ';'))
    gene_list = gsub(' ', '', gene_list) # Replacing any spaces with nothing
    
    interested.indices <- c()
    
    # Look for genes on the first row
    tmp <- intersect(genes_1, gene_list)
    for (j in 1 : length(tmp) ){
      interested.indices <- append(interested.indices, indices_1[tmp[j],1] : indices_1[tmp[j],2])
    }
    interested.indices <- ind_1$ix[interested.indices] # data.standard is already sorted by gene1 so it may not be necessary
    # unique(data_tmp_1[interested.indices,]$gene1)
    # unique(data_input[ind_1$ix[interested.indices],]$gene1)
    
    # Look for genes on the second row
    tmp <- intersect(genes_2, gene_list)
    interested.indices_tmp <- c()
    for (j in 1 : length(tmp) ){
      interested.indices_tmp <- append(interested.indices_tmp, indices_2[tmp[j],1] : indices_2[tmp[j],2])
    }
    # unique(data_tmp_2[interested.indices_tmp,]$gene2)
    # unique(data_input[ind_2$ix[interested.indices_tmp],]$gene2)
    
    interested.indices <- append(interested.indices, ind_2$ix[interested.indices_tmp])
    interested.indices <- unique(interested.indices) # data.standard[interested.indices,]
    
    ## Need to find only the positives that come from this complex (not this gene, as it can connect to other genes not part of the complex).
    # For negatives, it's fine.
    data_subset = entity.matrix[interested.indices,]   
    to.match = as.character(summary.standard$ID[i])
    
    # Now find the negative and positve part of the standard
    # neg.ind <- stri_isempty(data_subset$ID) 
    # identical(neg.ind, neg.ind.1) is TRUE!!! No need for 'stringi' package
    
    # All negatives are with ID = ''
    # "^\\s*$" asks for 0 or more (*) spaces (\\s) between beginning (^) and end ($) of string
    neg.ind <- grepl("^\\s*$", data_subset$ID)
    data.neg <- data_subset[neg.ind,]
    
    # Remove the pairs (ID) that don't contain the complex ID
    data.pos <- data_subset[!neg.ind, ]
    real.pos.ind <- sapply(sapply(strsplit(data.pos$ID, ';'), is.element, to.match), sum) == 1
    data.pos <- data.pos[real.pos.ind,]
    
    data_relevant <- rbind(data.neg, data.pos)
    # write.table(data_relevant[,c('true', 'predicted')], 'test_case_Prefoldin_33.txt', row.names = FALSE, sep = '\t')
    
    Perf.result = GenerateDataForPerfCurve(data_relevant$predicted, data_relevant$true, x.axis = 'recall', y.axis = 'precision')
    
    AUC.values[curr.ind] <- Perf.result$auc # Should be either zero or positive (check)
    diff.AUC[curr.ind] <- AUC.values[curr.ind] - Perf.result$y[length(Perf.result$y)] # differential AUPRC 
    
    setTxtProgressBar(pb, curr.ind/length(matched.ind)) # Progress bar update
    
    curr.ind <- curr.ind + 1
  }
  
  close(pb)
  
  # sort the data by their AUPRC values (highest to lowest)
  data.out <- data.frame(ID = summary.standard[matched.ind,]$ID, Name = summary.standard[matched.ind,]$Name, Length = summary.standard[matched.ind,]$Length, AUPRC = round(AUC.values,3), diff.AUPRC = diff.AUC, stringsAsFactors = FALSE)
  ind <- order (data.out$AUPRC, decreasing = TRUE) # why is ind$ix not giving the same data of ind$x!!
  data.out <- data.out[ind,]
  
  return (data.out)
}


#' Calculate the contribution of complex/PW etc. (num TP) at different precision cutoffs in a tabular format. This is different from Stepwise Contribution in that it doesn't find the smallest subset to explain all the TP, but takes all of them (and therefore some TPs are counted twice i.e. those are in multiple complexes)
#'
#' @param Pairs.in.data Identity of each entity ('ID'), 'Name', Genes insdie ('Gene'), 'Length'
#' @param list.of.precisions Precision cutoffs the user is interested in.
#' @param summary.standard Summary of the standard used: includes ID, Name, Genes, and Length
#'
#' @return output.complex.contribution -> Complex by precision matrix (where entries are number of TP contributed by the complex at that precision)
#' This is a entity ID * cutoff data.frame (with an extra column added for entity Name)
#' @export
#'  
GetContributionOfEntitiesAtPrecisions <- function (Pairs.in.data, list.of.precisions, summary.standard){
  
  ## *** Check input data format
  if((class(Pairs.in.data) != 'data.frame') & (class(summary.standard) != 'data.frame')){
    stop ("Check the input again ... it's supposed to be data.frame")
  }
  
  # Cacluate TP and Precision
  ind <- order(Pairs.in.data$predicted, decreasing = TRUE)
  Pairs.in.data <- Pairs.in.data[ind,]
  
  PR.values <- GenerateDataForPerfCurve(value.predicted = Pairs.in.data$predicted, value.true = Pairs.in.data$true, x.axis = 'TP', y.axis = 'precision')
  TP = PR.values$x
  Precision = PR.values$y
  ind.valid.precision <- (list.of.precisions >= min(Precision) & list.of.precisions <= max(Precision)) # Don't calculate for precisions out of range
  
  TP.count = cumsum(Pairs.in.data$true)
  
  ## Let's make a matrix of row = all complexes (IDs) and col = precision cutoffs
  # list.of.precisions <- c(0.0059, seq(0.1, 1, 0.025)) # Looks good
  Pos.Pairs.in.data <- Pairs.in.data[Pairs.in.data$true == 1,]
  unique.IDs <- unique(unlist(strsplit(paste(Pos.Pairs.in.data$ID, collapse = ';'), split = ';')))
  
  ID.cutoff.matrix <- matrix(0, nrow = length(unique.IDs), ncol = length(list.of.precisions)) # Initial TP contribution matrix
  rownames(ID.cutoff.matrix) <- as.numeric(unique.IDs) # unique.IDs
  
  for (i in seq(1, length(list.of.precisions), 1)) { # For all cutoffs provided
    # for (i in seq(1, 5, 1)) {
    # i = 20
    # print(paste('cutoff: ', list.of.precisions[i], sep = ': '))
    
    # If we don't have TP and Precision values for this, skip (will be zeros in output)
    if(!ind.valid.precision[i]){ next } 
    
    cutoff <- list.of.precisions[i]
    cand.ind <- which(Precision >= cutoff) # Last one (if > 1)
    tmp.ind <- which(TP.count == TP[cand.ind[length(cand.ind)]]) # TP[cand.ind[length(cand.ind)]] : number of TPs at this precision
    
    # Get the TP pairs and their association (to complex IDs)
    tmp.pairs <- Pairs.in.data[1 : tmp.ind[1], ] # Should we use the first one or the last one? I don't think it matters because, the number of TP is not increasing anyway here!s
    tmp.pairs <- tmp.pairs[tmp.pairs$true == 1, ]
    
    ## Find contribution of pairs for each entity (ID)
    for (j in seq(1, dim(tmp.pairs)[1], 1) ) {
      
      # each ID has value (as we are only keeping TPs)
      tmp.names <- unlist(strsplit(tmp.pairs$ID[j], ';'))
      
      # for each individual ID in the ID column
      for (k in tmp.names){
        ID.cutoff.matrix[k,i] <- ID.cutoff.matrix[k,i] + 1 # Initialized with 0
      }
    }
  }
  
  # Remove those IDs that didn't contribute to any TP pairs (in stepwise contribution)
  nonzero.cont.ind <- which(!apply(ID.cutoff.matrix, 1, sum) == 0)
  final.contribution.matrix <- ID.cutoff.matrix[nonzero.cont.ind, ]
  
  ## map the IDs to the names of the Complex and combine to a data.frame
  ia <- intersect(as.character(summary.standard$ID), row.names(final.contribution.matrix))
  # ind_name <- match(ia, as.character(summary.standard$ID)) # identical(ia, as.character(summary.standard[ind_name,]$ID))
  ind_std <- match(ia, as.character(summary.standard$ID))
  ind_contr <- match(ia, row.names(final.contribution.matrix))
  
  summary.contribution <- cbind(Name = summary.standard[ind_name,]$Name, as.data.frame(final.contribution.matrix[ind_contr,], stringsAsFactors = F))
  
  ## Add column names (if a colname starts with 0, r reads it as X0 !!)
  # names.prec <- gsub('.', '_', as.character(list.of.precisions), fixed = T) # fixed = T is important!!
  # colnames(summary.contribution)[2:dim(summary.contribution)[2]] <- names.prec
  # colnames(summary.contribution)[2:dim(summary.contribution)[2]] <- as.character(list.of.precisions)
  colnames(summary.contribution)[2:dim(summary.contribution)[2]] <- paste0('Precision_', as.character(cutoff.all))
  
  # Sort by the maximum Precision cutoff (where we have valid data)
  ind <- order(summary.contribution[, max(which(ind.valid.precision))], decreasing = T)
  summary.contribution <- summary.contribution[ind,]
  
  return (summary.contribution)
}



#' Concentrate the signal (TPs) into the minimum number of complexes that can explain the signal + get the contribution of the complexes at different precisions.
#'
#' @param Pairs.in.data Identity of each entity ('ID'), 'Name', Genes insdie ('Gene'), 'Length'
#' @param cutoff.all Precision cutoffs the user is interested in.
#' @param summary.standard Summary of the standard used: includes ID, Name, Genes, and Length
#'
#' @return output.stepwise.contribution -> Stepwise contribution matrix
#' This is a entity ID * cutoff data.frame (with an extra column added for entity Name)
#' @export
#'  
GetStepwiseContributionOfEntities <- function (Pairs.in.data, cutoff.all, summary.standard){
  
  ## *** Check input data format
  if((class(Pairs.in.data) != 'data.frame') & (class(summary.standard) != 'data.frame')){
    stop ("Check the input again ... it's supposed to be data.frame")
  }
  
  # Sort the data by score (prediction)
  ind <- order(Pairs.in.data$predicted, decreasing = TRUE)
  Pairs.in.data <- Pairs.in.data[ind,]
  
  # Just a little different as in matlab the first values are 0 (TP or PR.values$x) and Nan (Precision or PR.values$y)
  PR.values <- GenerateDataForPerfCurve(value.predicted = Pairs.in.data$predicted, value.true = Pairs.in.data$true, x.axis = 'TP', y.axis = 'precision')
  TP = PR.values$x
  Precision = PR.values$y
  TP.count = cumsum(Pairs.in.data$true)
  
  # *** Don't calculate for precisions out of range
  ind.valid.precision <- (cutoff.all >= min(Precision) & cutoff.all <= max(Precision))
  
  ## Let's make a matrix of row = all complexes (IDs) and col = precision cutoffs
  Pos.Pairs.in.data <- Pairs.in.data[Pairs.in.data$true == 1,]
  unique.IDs <- unique(unlist(strsplit(paste(Pos.Pairs.in.data$ID, collapse = ';'), split = ';')))
  
  # Initial TP contribution matrix
  ID.cutoff.matrix <- matrix(0, nrow = length(unique.IDs), ncol = length(cutoff.all)) 
  rownames(ID.cutoff.matrix) <- as.numeric(unique.IDs) # unique.IDs
  
  # Final TP contribution matrix
  final.contribution.matrix <- matrix(0, nrow = length(unique.IDs), ncol = length(cutoff.all)) # Stepwise TP contribution matrix
  rownames(final.contribution.matrix) <- unique.IDs
  
  for (i in seq(1, length(cutoff.all), 1)) {
  
    # If thsi precision cutoff is invalid (we don't have data for it), skip 
    # All contributions will remain zeros in output for these cutoffs 
    if(!ind.valid.precision[i]){ next } 
    
    print(paste('valid cutoff: ', cutoff.all[i], sep = ': '))
    
    cutoff = cutoff.all[i]
    cand.ind <- which(Precision >= cutoff)
    tmp.ind <- which(TP.count == TP[cand.ind[length(cand.ind)]]) # Looks good
    
    # Get the TP pairs and their association (to complex IDs)
    tmp.pairs <- Pairs.in.data[1 : tmp.ind[1], ]
    tmp.pairs <- tmp.pairs[tmp.pairs$true == 1, ] # Only keep the TPs
    
    ## Find contribution of pairs for each entity (ID)
    for (j in seq(1, dim(tmp.pairs)[1], 1) ) {
      
      # each ID has value (as we are only keeping TPs)
      tmp.names <- unlist(strsplit(tmp.pairs$ID[j], ';'))
      
      # for each individual ID in the ID column
      for (k in tmp.names){
        ID.cutoff.matrix[k,i] <- ID.cutoff.matrix[k,i] + 1 # Initialized with 0
      }
    }
    
    ## Find the maximum contribution at this cutoff and update the final.contribution.matrix
    max.ind <- which(ID.cutoff.matrix[,i] == max(ID.cutoff.matrix[,i], na.rm = TRUE))
    curr.id <- names(max.ind[1])
    final.contribution.matrix[curr.id,i] <- ID.cutoff.matrix[curr.id,i]
    
    # Until the TP pairs are all counted for, remove the ID that contributes the most, the next and so on
    tmp.pairs.small <- tmp.pairs
    id.ind <- 1
    store.sizes <- c(dim(tmp.pairs.small)[1])
    
    while (dim(tmp.pairs.small)[1] > 0) {
      
      ## Matching examples
      # The ID string will be present in the following form
      # \\b
      # \\b -> boundary 
      # ;? (one or zero ;) 
      # curr.id (actual ID string)
      # string = c("178", '123;178;124', '178;124', '123;178', "1178", "1788", '123;1178', '1178;124', '123;1178;124')
      # grepl("\\b;?178;?\\b", string, fixed = FALSE)
      # TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE
      curr.id.pattern <- paste0('\\b;?', curr.id, ';?\\b')
      log.ind <- unlist(lapply(tmp.pairs.small$ID, grepl, pattern = curr.id.pattern, fixed = FALSE))
      
      # Now remove pairs associated with that ID
      tmp.pairs.small <- tmp.pairs.small[!log.ind, ]
      # print(dim(tmp.pairs.small)[1])
      store.sizes <- c(store.sizes, dim(tmp.pairs.small)[1])
      
      if (dim(tmp.pairs.small)[1] == 0){
        break
      }
      
      unique.IDs.new <- unique(unlist(strsplit(paste(tmp.pairs.small$ID, collapse = ';'), split = ';')))
      ID.contribution <- matrix(0, nrow = length(unique.IDs.new), ncol = 1)
      names(ID.contribution) <- unique.IDs.new
      
      for (j in seq(1, dim(tmp.pairs.small)[1], 1) ) {
        
        # each ID has value (as we are only keeping TPs)
        tmp.names <- unlist(strsplit(tmp.pairs.small$ID[j], ';'))
        
        # for each individual ID in the ID column
        for (k in tmp.names){
          ID.contribution[k] <- ID.contribution[k] + 1 # Initialized with 0
        }
      }
      
      # Save ID and contribution
      max.ind <- which(ID.contribution == max(ID.contribution, na.rm = TRUE))
      curr.id <- names(max.ind[1])
      final.contribution.matrix[curr.id,i] <- ID.contribution[curr.id]
      
      id.ind <- id.ind + 1
    }
  
  }
  
  # Remove those IDs that didn't contribute to any TP pairs (in stepwise contribution)
  nonzero.cont.ind <- which(!apply(final.contribution.matrix, 1, sum) == 0)
  final.contribution.matrix <- final.contribution.matrix[nonzero.cont.ind, ]
  
  ## Now we have to map the IDs to the actual names of the entity (Complex)
  ia <- intersect(as.character(summary.standard$ID), row.names(final.contribution.matrix)) # Common complex IDs
  
  output.stepwise.contribution <- as.data.frame(final.contribution.matrix[ia, ])
  colnames(output.stepwise.contribution) <- paste0('Precision_', as.character(cutoff.all))
  
  # Add the name of complexes as the first column
  ind_name <- match(ia, as.character(summary.standard$ID)) # identical(ia, as.character(summary.standard[ind_name,]$ID))
  output.stepwise.contribution <- cbind(Name = summary.standard[ind_name,]$Name, output.stepwise.contribution, stringsAsFactors = F)
  # colnames(output.stepwise.contribution)[1] <- 'Name'
  
  # Sort by the maximum (valid) Precision cutoff
  ind <- order(output.stepwise.contribution[, max(which(ind.valid.precision))], decreasing = T)
  output.stepwise.contribution <- output.stepwise.contribution[ind,]
  
  # Convert factor to character
  # output.stepwise.contribution <- data.frame(lapply(output.stepwise.contribution, as.character), stringsAsFactors=FALSE)
  
  return (output.stepwise.contribution)
}



#' Remove a complex (or complexes) to re-evaluate the data. This removes all the positive pairs associated with them complex (thereby effectively removing the complex).
#' @param data.standard input data: either co-annotation standard or output of CalculatePredictionAndTrueOnLibraryProfiles as a data.frame (need the 'ID' column for both)
#' @param ids the complex ids to remove
#' @param replace
#'  Way1 (replace = false): Remove the positive pairs associated to ids.
#'  This reduces the size of the data a bit.
#'  Way2 (replace = true): Convert the positive pairs to negatives. 
#'  This will maintain the size (number of pairs) of the data.
#' @export
#' 
getSubsetOfCoAnnRemoveIDs <- function(data.standard, ids, replace = FALSE){
  
  ## *** Check input data format
  if(class(data.standard) != 'data.frame'){
    stop ("data.standard is supposed to be data.frame")
  }
  if(class(ids) != 'character'){
    stop ("ids is supposed to be a character vector")
  }
  
  ## Get the positive examples (TPs)
  if (sum(grepl('true', colnames(data.standard)))){
    log.ind.pos = which(data.standard$true == 1)  
  } else if (sum(grepl('is_annotated', colnames(data.standard)))){
    log.ind.pos = which(data.standard$is_annotated == 1)  
  } else{
    stop ("ill-formated input: no column named 'is_annotated' or 'true' ")
    return (NULL)
  }
  data.input = data.standard[log.ind.pos,] # 39583
  
  interested.indices = c()
  for (complex.id.str in ids){ # '320'
    # complex.id.str <- ids
    
    # Find those that contain '320'
    log.ind.1 = which(grepl(complex.id.str, data.input$ID) == TRUE) # 3021
    # data.input[log.ind.1, ]
    
    # Among those from log.ind.1, there maybe some that contain '320'
    # as subset ('6320' or '3204' for example). We want to keep those.
    
    # Add back the data that contains '6320' but not '320'
    # IDs <- lapply(strsplit(data.input$ID[log.ind.1], ';'), %in%, complex.id.str)
    id.membership <- lapply(strsplit(data.input$ID[log.ind.1], ';'), is.element, complex.id.str)
    log.ind.2 <- which(unlist(lapply(id.membership, sum)) == 1)
    interested.indices = c(interested.indices, log.ind.1[log.ind.2])
    
    data.input[interested.indices, ]
  }
  
  interested.indices = unique(interested.indices) # 3003
  length(interested.indices)

  data.standard[log.ind.pos[interested.indices], ]
  
  # Final data
  if (replace == FALSE){ # Way1: Remove the positive examples associated with '320'
    data.output = data.standard
    data.output <- data.output[-log.ind.pos[interested.indices], ]
  } 
  else { # Way 2: Convert the positive (1) examples associated with '320' to negative (0)
    data.output = data.standard
    if (sum(grepl('true', colnames(data.standard)))){
      data.output$true[log.ind.pos[interested.indices]] = 0
    } else{
      data.output$is_annotated[log.ind.pos[interested.indices]] = 0
    }
  }
 
  return (data.output) 
}



#' Remove a set of genes to re-evaluate the data. Should be used only when we can't use the sister function getSubsetOfCoAnnRemoveIDs or getSubsetOfCoAnnRemovePairs.
#' @param data_standard A co-annotation standard (data.frame)
#' @param data_subset Output of CalculatePredictionAndTrueOnLibraryProfiles with index column (data.frame)
#' @param gene_list a vector of genes to remove. This will, for positive examples, remove all the pairs that include the gene. (character vector)
#' @param replace
#'  Way1 (replace = false): Remove the positive pairs associated to ids.
#'  This reduces the size of the data a bit.
#'  Way2 (replace = true): Convert the positive pairs to negatives. 
#'  This will maintain the size (number of pairs) of the data.
#' @export
#' 
getSubsetOfCoAnnRemoveGenes <- function(data_standard, data_subset, gene_list, replace = FALSE){
  
  ## *** Check input data format
  if( (class(data_standard) != 'data.frame') | (class(data_subset) != 'data.frame')){
    stop ("data_standard is supposed to be data.frame")
  }
  if(class(gene_list) != 'character'){
    stop ("ids is supposed to be a character vector")
  }
  
  ## *** Check if index column is given for data_subset
  if (!sum(grepl('index', colnames(data_subset)))){
    stop('No indices given for data_subset .. we need this!!')
  }
  
  ## Get the positive examples (TPs)
  if (sum(grepl('true', colnames(data_subset)))){
    ind_pos = which(data_subset$true == 1)  
  } else if (sum(grepl('is_annotated', colnames(data_subset)))){
    ind_pos = which(data_subset$is_annotated == 1)  
  } else{
    stop ("ill-formated input: no column named 'is_annotated' or 'true' ")
    return (NULL)
  }
  
  tmp = data_subset[ind_pos,]
  # identical(data_standard$ID[tmp$index], tmp$ID) # Should be true: Sanity check
  # Now get the gene symbols from data_standard and merge
  row.names(tmp) <- tmp$index # The indices are the row identifiers in data_standard
  data.input <- cbind(data_standard[tmp$index,c('gene1', 'gene2')], tmp)
  
  ## Sort the data by first gene (and group them)
  ind <- order(data.input$gene1)
  data.tmp <- data.input[ind,]
  gene.indices <- GroupUniqueElements(data.tmp$gene1)
  
  genes_overlapping <- sort(intersect(gene_list, row.names(gene.indices)))
  interested_indices = c()
  for (i in genes_overlapping){
    interested_indices = c(interested_indices, gene.indices[i,1] : gene.indices[i,2])
  }
  interested_indices = ind[interested_indices]
  
  ## Sort the data by second gene (This part is usually unnecessary as our data is sorted and there are no duplicate pairs)
  ind <- order(data.input$gene2)
  data.tmp <- data.input[ind,]
  gene.indices <- GroupUniqueElements(data.tmp$gene2)
  
  genes_overlapping <- sort(intersect(gene_list, row.names(gene.indices)))
  interested_indices_tmp = c()
  for (i in genes_overlapping){
    interested_indices_tmp = c(interested_indices_tmp, gene.indices[i,1] : gene.indices[i,2])
  }
  
  interested_indices = c(interested_indices, ind[interested_indices_tmp])
  interested_indices = unique(interested_indices) # 3003
  length(interested_indices)
  
  # data.input[interested_indices,]
  # data_subset[ind_pos[interested_indices], ]
  ## *** The problem of this is it can remove some pairs where one gene come from the complex, and the other from somewhere else!!!
  if (FALSE){
    discrepancy.ind <- setdiff(data_subset[ind_pos[interested_indices], 'index'], data.standard[log.ind.pos[interested.indices], 'index'])
    tmp.match <- intersect(data_subset$index, discrepancy.ind)
    tmp.ind <- match(tmp.match, data_subset$index)
    data_standard[data_subset[tmp.ind,'index'], ]  
  }
  
  # Final data
  if (replace == FALSE){ # Way1 (default): Remove the positive examples associated with '320'
    data.output <- data_subset
    data.output <- data.output[-ind_pos[interested_indices], ]
  } 
  else { # Way 2: Convert the positive examples associated with '320' to negative
    data.output <- data_subset
    
    if (sum(grepl('true', colnames(data.standard)))){
      data.output$true[ind_pos[interested_indices]] = 0
    } else{
      data.output$is_annotated[ind_pos[interested_indices]] = 0
    }
  }
  
  return (data.output) 
}



#' Remove gene-pairs to re-evaluate the data. This is probably more accurate (albeit time consuming) than getSubsetOfCoAnnRemoveGenes
#' @param data_standard A co-annotation standard (data.frame)
#' @param data_subset Output of CalculatePredictionAndTrueOnLibraryProfiles with index column (data.frame)
#' @param gene_list Any co-functionality between genes in this list will be removed.
#' @param replace
#'  Way1 (replace = false): Remove the positive pairs associated to ids.
#'  This reduces the size of the data a bit.
#'  Way2 (replace = true): Convert the positive pairs to negatives. 
#'  This will maintain the size (number of pairs) of the data.
#' @export
#' 
getSubsetOfCoAnnRemovePairs <- function(data_standard, data_subset, gene_list, replace = FALSE){
  
  ## *** Check input data format
  if(class(data_subset) != 'data.frame'){
    stop ("data_subset is supposed to be data.frame")
  }
  if(class(gene_list) != 'list'){
    stop ("gene_list is supposed to be a list of genes")
  }
  
  if(class(data_standard) == 'list'){
    ## Get the groupings
    gene.indices <- data_standard$gene.indices # sorted by entrez ID
    unique.genes.std.entrez <- row.names(gene.indices)
    
    ## Get the mappings
    mapping <- data_standard$mapping
    unique.genes.std.symbol <- mapping[row.names(gene.indices)]
    
    # Use the mappings to change gene_list
    mapped_list <- gene_list
    for (i in 1:length(mapped_list) ){
      common_genes <- intersect(mapped_list[[i]], mapping)
      ind_genes_in_mapping <- match(common_genes, mapping)
      mapped_list[[i]] <- as.integer(names(mapping)[ind_genes_in_mapping])
    }
    # identical(unname(mapping[as.character(mapped_list[[1]])]), gene_list[[1]]) # TRUE
    
    ## Main data
    data.relevant <- data_standard$data
    
  } else if(class(data_standard) == 'data.frame'){
    data.relevant <- data_standard
    mapped_list <- gene_list
    
  } else{
    stop('data_standard is supposed to be a data.frame or a list')
  }
  
  ## *** Check if gene pairs are provided for data_standard
  pair_sum <- sum(grepl('gene1', colnames(data.relevant))) + 
    sum(grepl('gene2', colnames(data.relevant)))
  if (pair_sum != 2){
    stop('data_standard: gene pair names are not provided !!')
  }
  
  ## ===================================================
  ## Get the positive pairs from the provided inputs
  if (!is.null(data_subset)){ # data_subset is provided
    if (!sum(grepl('index', colnames(data_subset)))){
      stop('No indices given for data_subset .. we need this!!')
    }
    if (sum(grepl('true', colnames(data_subset)))){
      log_ind_pos = which(data_subset$true == 1) # log_ind_pos: Need later ***
      tmp = data_subset[log_ind_pos,]
      row.names(tmp) <- tmp$index # The indices are the row identifiers in data_standard
      
      data_input <- cbind(data.relevant[tmp$index,c('gene1', 'gene2')], tmp)
    } else{
      stop ("ill-formated input for data_subset: no column named 'true' ")
      return (NULL)
    }
    
    data.output <- data_subset # Initializing output
    
  } else{ # data_subset is NULL / not provided
    if (sum(grepl('is_annotated', colnames(data.relevant)))){
      log_ind_pos = which(data.relevant$is_annotated == 1)
      data_input <- data.relevant[log_ind_pos, ]
      
      data.output <- data.relevant # Initializing output
    } else{
      stop ("ill-formated input for data_standard: no column named 'is_annotated' ")
      return (NULL)
    }
  }

  ## Sort the data by first gene and group
  ind_saved <- order(data_input$gene1) # ind_saved: Need later ***
  data_tmp <- data_input[ind_saved,]
  
  ind_source <- GroupUniqueElements(data_tmp$gene1)
  genes_source <- row.names(ind_source)
  
  
  ## -------- Construct the list of unique pairs to remove / replace --------
  gene_first <- c()
  gene_second <- c()
  count <- 0
  
  for (i in 1:length(mapped_list) ){
    genes_in_complex <- sort(mapped_list[[i]])
    count <- count + choose(length(genes_in_complex), 2)
    
    for (j in 1: (length(genes_in_complex)-1) ){ # The -1 is necessary
      for (k in (j+1):length(genes_in_complex) ){
        gene_first = c(gene_first, genes_in_complex[j])
        gene_second = c(gene_second, genes_in_complex[k])
      }  
    }
  }
  
  # Only keep the unique pairs
  gene_combined <- paste(gene_first, gene_second, sep = '_')
  unique.index <- which(!duplicated(gene_combined))
  pair_list <- as.data.frame(cbind(gene_first[unique.index], gene_second[unique.index]), stringsAsFactors = FALSE)
  colnames(pair_list) <- c('first', 'second') # 969 candidate pairs for removal (ETCI) - good
  
  # Now sort the pairs and group them
  ind <- order(pair_list$second) # by second gene
  pair_list <- pair_list[ind,]
  ind <- order(pair_list$first) # then by first gene
  pair_list <- pair_list[ind,]
  
  ind_cand <- GroupUniqueElements(pair_list$first) # group by first gene
  genes_cand <- row.names(ind_cand) # ETCI (45) - good
  
  ## *** Now find the pairs (indices) to remove in the source/annotation data
  genes_common <- intersect(genes_source, genes_cand) # 70
  ia <- match(genes_common, genes_source)
  ic <- match(genes_common, genes_cand)
  
  interested_indices <- c() # Good until here?
  for (i in seq_along(ia)){ # i in 1 : length(ia) will fail if ia is empty
    sec_genes_cand <- pair_list$second[ind_cand[ic[i],1]:ind_cand[ic[i],2]]
    
    src_indices <- ind_source[ia[i],1] : ind_source[ia[i],2]
    sec_genes_source <- data_tmp$gene2[src_indices]
    
    genes_tmp <- intersect(sec_genes_source, sec_genes_cand)
    ia_source <- match(genes_tmp, sec_genes_source)
    
    interested_indices = c(interested_indices, src_indices[ia_source]) # *** And this one
  }
  (length(interested_indices))

  # data_tmp[interested_indices,] # unique(data_tmp[interested_indices,1])
  # one of the following are going to be true
  # identical(data_subset[log_ind_pos[ind_saved[ interested_indices]],'is_annotated'], data_tmp[interested_indices,'is_annotated'])
  # identical(data_subset[log_ind_pos[ind_saved[ interested_indices]],'true'], data_tmp[interested_indices,'true'])
  
  # Check for unexpected output
  if (length(interested_indices) == 0){
    return (data.output)
  }
  
  # mapping[as.character(unique(data.relevant[data.output[log_ind_pos[ind_saved[interested_indices]], 3],1]))]
  
  # Final data
  if (replace == FALSE){ # Way1 (default): Remove the positive examples associated with '320'
    data.output <- data.output[-log_ind_pos[ind_saved[interested_indices]], ]
  } else { # Way 2: Convert the positive examples associated with '320' to negative
    if (sum(grepl('true', colnames(data.relevant)))){
      data.output$true[log_ind_pos[ind_saved[interested_indices]]] <- 0
    } else{
      data.output$is_annotated[log_ind_pos[ind_saved[interested_indices]]] <- 0
    }
  }
  
  return (data.output) 
}
