## Generate different plots for visualization
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


#' Plot Performance Curves on Profile Similarities
#'
#' @param pred.ca a list of elements each contain two lists predicted (score) and true (co-anotation)
#' @param subsample Subsamples from the PR space to reduce the size of the plotted data. 
#' @param neg.to.pos assign TRUE to sort from neg to positive (pos -> neg is the default)
#' @param type.plot can be either 'log' - semilog (x axis) or 'normal' plot
#' @param fig.title title of the figure (default nothing)
#' @param fig.labs x and y axis labels - default: c('TP', 'Precision')
#' @param legend.names names of legends
#' @param legend.color colors of legends
#' @param save.figure set to TRUE to save the figure in a pdf file (fig.title is used as the name)
#' @param outfile.name the name of the output file(figure)
#' @param outfile.type type of figure to save - 'pdf' (default) or 'png'
#' @param is.bgdline set to TRUE if we want to plot a background (reference) line for precision (y-axis). If plotting multiple plots, setting this to false (default) probably makes more sense unless all of those share the same background.
#'
#' @return
#' 
#' @examples
#' 
#' @export

PlotPRSimilarity <- function(pred.ca, subsample = FALSE,
                             neg.to.pos = FALSE, 
                             type.plot = 'log', is.bgdline = FALSE,
                             provided.xlim = NULL, provided.ylim = NULL,
                             legend.names = NULL, legend.color = c('blue'),
                             legend.ltype = NULL, box.type = 'L',
                             fig.title = NULL, fig.labs = c('TP', 'Precision'),
                             save.figure = FALSE, 
                             outfile.name = 'test_PR_sim', outfile.type = 'pdf') {
  
  ## *** Check input data format
  if(class(pred.ca) != 'list'){
    stop ('A list of list is expected as input ... ')
  } else{ # check if we have the correct data inside the list
    
    corr_columns <- sum(grepl('true', names(pred.ca[[1]]))) + sum(grepl('predicted', names(pred.ca[[1]])))
    if(corr_columns != 2){
      stop (" All list element should contain a 'true' and a 'predicted' list ... ")
    }
  }
  
  if (length(legend.color) < length(pred.ca)){
    warning('Color not provided for each curve !! Making our own ')
    if (length(pred.ca) == 1){
      legend.color = c('blue')
    }
    else{
      legend.color = palette(rainbow(n = length(pred.ca)))
    }
  }
  
  
  ## *** Make the data to plot on the x (TP) and y (Precision) axis
  for (i in 1 : length(pred.ca)){
    
    test <- GenerateDataForPerfCurve(value.predicted = pred.ca[[i]]$predicted, value.true = pred.ca[[i]]$true, neg.to.pos = neg.to.pos, x.axis = 'TP', y.axis = 'precision')
    
    ## ** If we want a faster plot (less points)
    if (subsample){
      
      # Take the last precision value for the same TP
      unique.index <- which(!duplicated(test$x))
      unique.index <- c(unique.index - 1, length(test$x))
      unique.index <- unique.index[-1]
      
      small_x <- test$x[unique.index]
      small_y <- test$y[unique.index]
      
      # Now subsample from the whole space (top 100, and then 100 from each 10^ range)
      highest.power <- ceiling(log10(length(small_x)))
      
      if (highest.power > 3){
        # All from 1 to 100
        keep.index <- c(1:100)
        
        # 100 from each interval (i.e 1001 to 10000)
        for (pow in 3 : (highest.power - 1)){
          keep.index <- c(keep.index, c(round(seq(10^(pow-1)+1, 10^pow, length = 100))))
        }
        max.available <- min(100, length(small_x) - (10^pow))
        keep.index <- c(keep.index, 
                        c(round(seq(10^pow+1, length(small_x), 
                                    length = max.available))))
        
        test$x <- small_x[keep.index]
        test$y <- small_y[keep.index]
        
      } else{ # If the number of points are < 1000, no need for subsampling
        test$x <- small_x
        test$y <- small_y
      }
    }
    
    if (i == 1){
      plot.data <- list(list(x = test$x, y = test$y))
    } else{
      plot.data <- append(plot.data, list(list(x = test$x, y = test$y)))
    }
  } # end for
  
  
  ## *** Calculate the max y-lim and x-lim for the plots (if not provided)
  #  In R, we have to this beforehand as we are going to plot multiple lines
  plot.xlim = -1
  plot.ylim = -1
  
  for (i in 1 : length(plot.data)) {
    tmp <- plot.data[[i]]
    
    ind.10 <- which(tmp$x == 9) # Remove data for TP < 10
    
    if (length(ind.10) > 0){ 
      # The first 10 points are excluded from plotting
      x <- max(tmp$x[-(1:ind.10[length(ind.10)])])
      y <- max(tmp$y[-(1:ind.10[length(ind.10)])])
    } else{ # What if it's not there?
      x <- max(tmp$x)
      y <- max(tmp$y)
    }
    
    if (x > plot.xlim){
      plot.xlim <- x
    }
    if (y > plot.ylim){
      plot.ylim <- y
    }
  }
  
  if(!is.null(provided.xlim)) {plot.xlim <- provided.xlim}
  if(!is.null(provided.ylim)) {plot.ylim <- provided.ylim}
  
  if (type.plot == 'log'){
    plot.xlim <- 10 ^ ceiling(log10(plot.xlim))  
  }
  
  plot.ylim <- round(plot.ylim + 0.01, 2)
  
  ## *** Save to an output file
  if (save.figure == TRUE){
    if(outfile.type == 'png'){
      png(paste(outfile.name, ".png", sep = ''), width = 7, height = 5, units="in", res = 300)
    }else{
      pdf(paste(outfile.name, ".pdf", sep = '') )   
    }
  }
  
  for (i in 1 : length(plot.data)) {
    # print(i)
    
    tmp <- plot.data[[i]]
    ind.10 <- which(tmp$x == 9) # Remove data for TP < 10
    
    if (length(ind.10) > 0){ # What if it's not there?
      # The first 10 points are excluded from plotting
      data.x.axis <- tmp$x[-(1:ind.10[length(ind.10)])]
      data.y.axis <- tmp$y[-(1:ind.10[length(ind.10)])]
    } else{
      data.x.axis <- tmp$x
      data.y.axis <- tmp$y
    }
    
    if (i == 1){ # For the first time
      if (type.plot == 'log'){
        
        plot(data.x.axis, data.y.axis, 
             xlim = c(10, plot.xlim), ylim = c(0, plot.ylim),
             log = 'x', type = 'l', main = fig.title,
             bty = box.type, # 'n' -> no box, nothing - all boxes
             xlab = fig.labs[1], ylab = fig.labs[2], 
             lwd = 2, col = legend.color[i], lty = legend.ltype[i],
             cex.lab = 1.4, cex.main = 1.4, cex.axis = 1.4,
             font.main = 1, xaxt="n") # Not plotting x axis tick labels
        
        # Adding customized xtick labels
        # x.axis.ticks <- seq(1, (floor(log10(plot.xlim))), 1)
        x.axis.ticks <- seq(1, log10(plot.xlim), 1)
        x.tick.labels <- as.expression(lapply(x.axis.ticks, function(E) bquote(10^.(E))))
        axis(1, at=10^(x.axis.ticks), 
             labels=as.expression(lapply(x.axis.ticks, function(E) bquote(10^.(E)))), 
             cex.axis = 1.4)
      }
      else{
        # Normal plot
        plot(data.x.axis, data.y.axis, xlim = c(10, plot.xlim), ylim = c(0, plot.ylim), 
             type = 'l', main = fig.title,
             bty = box.type, # 'n' -> no box, nothing - all boxes
             xlab = fig.labs[1], ylab = fig.labs[2], font.main = 1,
             lwd = 2, col = legend.color[i], lty = legend.ltype[i],
             cex.lab = 1.4, cex.main = 1.2, cex.axis = 1.4)
      }
    }
    else{ # For every other time
      lines(data.x.axis, data.y.axis, col=legend.color[i], lwd = 2, lty = legend.ltype[i]) 
    }
  }
  
  if (!is.null(legend.names)){
    legend("topright", legend = legend.names, 
           fill = legend.color, lwd = 2, lty = legend.ltype[i],
           bty = box.type, # to remove the box
           cex = 1.2, text.col = "black", horiz = F)
    
    # Printing some warning, just in case!
    if (length(pred.ca) > length(legend.names)){
      warning('Legend not provided for all curves !!')
    }
    if (length(legend.color) != length(legend.names)){
      warning('Number of legends and colors are not equal !!!')
    }
  }
  
  if (is.bgdline){
    abline(h = min(data.y.axis), col = 'black', lty=2)
  }
  
  if (save.figure == TRUE){
    dev.off()
  }
}


#' Plot performance curves on direct interactions (Negative vs Positive)
#'
#' @param plot.data a list of two lists names x and y (gives the data on x and y axis)
#' @param type.plot can be either semilog (on x axis) (default) or regular
#' @param fig.title title of the figure
#' @param fig.labs x and y axis labels
#' @param outfile.name the name of the output file(figure)
#' @param outfile.type type of figure to save - 'pdf' (default) or 'png'
#' @param save.figure if TRUE, saves the figure after this run (uses the same name as fig.title)
#'
#' @examples
#' 
#' @export
PlotPRDirect <- function(plot.data, type.plot = 'log', 
                         fig.title = NULL, fig.labs = c('TP', 'Precision'),
                         outfile.name = 'test_DI', outfile.type = 'pdf', 
                         save.figure = FALSE) {
  
  ## *** Calculate the PR data to plot for positive and negative interactions
  score.pos.int <- plot.data$data$Score [plot.data$data$Score > 0]
  true.pos.int <- plot.data$data$True  [plot.data$data$Score > 0]
  pos.PR <- GenerateDataForPerfCurve(value.predicted = score.pos.int, value.true = true.pos.int, x.axis = 'TP', y.axis = 'precision', neg.to.pos = FALSE)
  
  score.neg.int <- plot.data$data$Score [plot.data$data$Score < 0]
  true.neg.int <- plot.data$data$True [plot.data$data$Score < 0]
  neg.PR <- GenerateDataForPerfCurve(value.predicted = score.neg.int, value.true = true.neg.int, x.axis = 'TP', y.axis = 'precision', neg.to.pos = TRUE)
  
  # pred.ca <- list(positive = list(true = true.pos.int, predicted = score.pos.int))
  
  ## Calculate the max y-lim and x-lim for the plots 
  ind.pos.10 <- which(pos.PR$x == 9)
  ind.neg.10 <- which(neg.PR$x == 9)
  
  if (length(ind.pos.10) > 0){ # What if it's not there?
    # Points until TP of 10 are excluded
    plot.xlim = max(max(pos.PR$x[-(1:ind.pos.10[length(ind.pos.10)])]), 
                    max(neg.PR$x[-(1:ind.neg.10[length(ind.neg.10)])]))
  } else{
    plot.xlim = max(max(pos.PR$x), max(neg.PR$x))
  }
  
  if (length(ind.neg.10) > 0){
    plot.ylim = max(max(pos.PR$y[-(1:ind.pos.10[length(ind.pos.10)])]), 
                    max(neg.PR$y[-(1:ind.neg.10[length(ind.neg.10)])]))
  } else{
    plot.ylim = max(max(pos.PR$y), max(neg.PR$y))
  }
  
  plot.xlim <- 10 ^ ceiling(log10(plot.xlim)) # Getting the next log10 boundary
  plot.ylim <- round(plot.ylim + 0.01, 2)
  
  
  ## Make arrangements to save the plot to an output file
  colors.direct <- c('#FFFF00', '#00A4FF') # (pos, neg)
  if (save.figure == TRUE){
    if (outfile.type == 'png'){
      png(paste(outfile.name, ".png", sep = ''), width = 4, height = 4, units="in", res = 300) 
    }else{
      pdf(paste(outfile.name, ".pdf", sep = '') )   
    }
  }
  
  ## ------------- 1. Plotting direct positives -------------
  if (length(ind.pos.10) > 0){
    data.x.axis <- pos.PR$x[-(1:ind.pos.10[length(ind.pos.10)])]
    data.y.axis <- pos.PR$y[-(1:ind.pos.10[length(ind.pos.10)])]
  } else{
    data.x.axis <- pos.PR$x
    data.y.axis <- pos.PR$y
  }
  
  # Adding customized xtick labels
  x.axis.ticks <- seq(1, (floor(log10(plot.xlim))), 1)
  x.tick.labels <- as.expression(lapply(x.axis.ticks, function(E) bquote(10^.(E))))
  
  if (type.plot == 'log'){
    plot(data.x.axis, data.y.axis, 
         xlim = c(10, plot.xlim), ylim = c(0, plot.ylim), 
         log = 'x', type = 'l', bty = 'L', main = fig.title,
         xlab = fig.labs[1], ylab = fig.labs[2], lwd = 2, 
         cex.lab = 1.2, cex.main = 1.2, cex.axis = 1.2,
         col = colors.direct[1], xaxt='n') 
    axis(1, at=10^(x.axis.ticks), 
         labels=as.expression(lapply(x.axis.ticks, function(E) bquote(10^.(E)))),
         cex = 1.2)
  } else{
    # Normal plot
    plot(data.x.axis, data.y.axis, xlim = c(10, plot.xlim), ylim = c(0, plot.ylim), 
         type = 'l', 
         bty = "L", # 'n' -> no box, nothing - all boxes
         xlab = fig.labs[1], ylab = fig.labs[2], lwd = 2, col = colors.direct[1], cex = 1.2)
  }
  
  ## ------------- 2. Overlaying direct negatives -------------
  if (length(ind.neg.10) > 0){
    data.x.axis <- neg.PR$x[-(1:ind.neg.10[length(ind.neg.10)])]
    data.y.axis <- neg.PR$y[-(1:ind.neg.10[length(ind.neg.10)])]
  } else{
    data.x.axis <- neg.PR$x
    data.y.axis <- neg.PR$y
  }
  lines(data.x.axis, data.y.axis, col=colors.direct[2], lwd = 2) 
  
  
  legend("topright", legend = c('Positive', 'Negative'), fill = colors.direct, 
         cex = 1, bty = "n", # remove box
         text.col = "black", horiz = F)
  
  # When all the plotting is done (Need to work on this)
  if (save.figure == TRUE){
    # par (new = F)
    dev.off()
  }
}


#' Plot a scatter plot of contribution of all complexes
#'
#' @param plot.data a list or a data.frame with columns 'Name', 'Length', and 'AUPRC'
#' @param fig.title title of the figure
#' @param fig.labs labels for x and y axis
#' @param show.text set TRUE to show the names of plotted complexes
#' @param show.cutoffs set TRUE to show the AUPRC and size cutoff chosen
#' @param save.figure: if TRUE, saves the figure as a pdf (using fig.title as name)
#' @param outfile.type type of figure to save - 'pdf' (default) or 'png'
#' @param outfile.name the name of the output file(figure)
#' 
#' @export
#'
#' @examples

PlotContributionScatter <- function(plot.data, 
                                    length.cutoff = 30, AUPRC.cutoff = 0.4, 
                                    fig.title = NULL, fig.labs = c('AUPRC', 'Size'),
                                    show.text = FALSE, show.cutoffs = FALSE,
                                    save.figure = FALSE, 
                                    outfile.type = 'pdf', outfile.name = 'test_scatter') {
  
  ## *** Check if we have the right columns in plot.data
  corr_columns <- sum(grepl('Name', names(plot.data))) + sum(grepl('Length', names(plot.data))) + sum(grepl('AUPRC', names(plot.data)))
  if(corr_columns < 3){
    stop ("plot.data should contain 'Name', 'Length', and 'AUPRC' values ... ")
  }
  
  ## *** Group by complex size and AUPRC values
  #  TODO: What if one of these are empty?
  ind_hi_low <- which((plot.data$Length > length.cutoff) & (plot.data$AUPRC <= AUPRC.cutoff))
  ind_low_hi <- which((plot.data$Length <= length.cutoff) & (plot.data$AUPRC > AUPRC.cutoff))
  ind_hi_hi  <- which((plot.data$Length > length.cutoff) & (plot.data$AUPRC > AUPRC.cutoff))
  
  name_size_hi_auprc_low <- plot.data$Name[ind_hi_low] # blue
  name_size_low_auprc_hi <- plot.data$Name[ind_low_hi] # green
  name_size_high_auprc_hi <- plot.data$Name[ind_hi_hi] # purple
  
  ## *** Assign colors
  abm <- list(name_size_hi_auprc_low, name_size_low_auprc_hi, name_size_high_auprc_hi)
  pcol <- rep("white", length(plot.data$AUPRC))
  
  for(i in 1:length(abm)) {
    pcol[plot.data$Name %in% abm[[i]]] <- c("#6baed6","#74c476", "#630098")[i]
  }
  
  if (save.figure == TRUE){
    if(outfile.type == 'png'){
      png(paste(outfile.name, ".png", sep = ''), width = 4, height = 4, units="in", res = 300)
    }else{
      pdf(paste(outfile.name, ".pdf", sep = '')) 
      # pdf(file = "contribution_scatter.pdf", width = 3.4, height = 3.8, useDingbats = F)
    }
  }
  
  y.max <- max(plot.data$Length, na.rm = T) + 10
  
  # Scatter plot
  plot(plot.data$AUPRC, plot.data$Length, las = 1, 
       # xlim = c(0, 1), ylim = c(0, 150),
       xlim = c(0, 1), ylim = c(0, y.max),
       xlab = fig.labs[1], ylab = fig.labs[2], 
       pch = 21, bg = pcol, bty = "n", lwd=.33, cex = 1.2, main = fig.title)
  
  # Adding lines to show cutoffs
  if (show.cutoffs){
    abline(h = length.cutoff, col = 'gray60', lty=2)
    # text(0.8, length.cutoff, paste0("Size = ", length.cutoff), col = "gray60", adj = c(0, -.1))
    abline(v = AUPRC.cutoff, col = 'gray60', lty=2)
    # text(AUPRC.cutoff, 120, paste0("AUPRC = ", AUPRC.cutoff), col = "gray60", adj = c(0, -.1), srt = 90)  
  }
  
  # Add txt
  if (show.text){
    if (length(ind_hi_low) > 0) text(plot.data$AUPRC[ind_hi_low], plot.data$Length[ind_hi_low], labels=plot.data$Name[ind_hi_low], cex=0.6)
    if (length(ind_low_hi) > 0) text(plot.data$AUPRC[ind_low_hi], plot.data$Length[ind_low_hi], labels=plot.data$Name[ind_low_hi], cex=0.6)
    if (length(ind_hi_hi) > 0)text(plot.data$AUPRC[ind_hi_hi], plot.data$Length[ind_hi_hi], labels=plot.data$Name[ind_hi_hi], cex=0.6) 
  }
  
  # Add legends
  legend("topright", 
         legend = c(expression('Size'[hi]*', AUPRC'[lo]), 
                    expression('Size'[lo]*', AUPRC'[hi]),
                    expression('Size'[hi]*', AUPRC'[hi]) ), 
         col = c("#6baed6","#74c476", "#630098"), pch = 19, # to show solid circles
         cex = 1, text.col = "black", horiz = F)
  
  # plot(1:10, xlab=expression('Size'[hi]*', AUPRC'[low]))
  
  if (save.figure == TRUE){
    dev.off()
  }
}



#' Plot contribution structure (diversity) of the complexes (a muller plot)
#'
#' @param plot.data a complex (unique) vs precision matrix where each element denote number of TP at that combination
#' @param cutoff.all all the precision cutoffs used here (the same as used for plot.data) 
#' @param min.pairs all the precision cutoffs used here (the same as used for plot.data) 
#' @param num.complex.to.show How many complexes we want to show (everything else will be put to others)?
#' @param list.of.complexes.to.show Use this list of complex to show on the plot, regardless of their ranking.
#' @param alternative.names Provide a list of names (corresponds to list.of.complexes.to.show) to use as alternative to original names
#' @param min.precision.cutoff How far down should we go in precision cutoff to calcualte contribution of complexes? Default is 0.5, meaning we calculate contributions starting from the highest precision (where we have at least min.pairs) to the min.precision.cutoff and take a mean contribution to rank the complexes.
#' @param ccol colors for the top complexes highlighted (top 10 contributing complexes are colored by default)
#' @param show.legend default TRUE. Set to FALSE if we don't want to print legends.
#' @param fig.title title in case we want to save the image
#' @param fig.labs x and y axis labels
#' @param save.figure do we want to save or just plot
#' @param outfile.name name of the output file
#' @param outfile.type type of figure to save - 'pdf' (default) or 'png'
#' 
#' @export
#'
#' @examples
#

PlotContributionStructure <- function(plot.data, cutoff.all = NULL, 
                                      min.pairs = 10, min.precision.cutoff = 0.5, 
                                      num.complex.to.show = 10, show.legend = TRUE,
                                      list.of.complexes.to.show = NULL, 
                                      alternative.names = NULL,
                                      ccol = NULL, y.lim = NULL, fig.title = NULL, 
                                      fig.labs = c('Fraction of TP', 'Precision'), 
                                      outfile.name = 'test_cont_str', outfile.type = 'pdf',
                                      save.figure = FALSE) {
  
  ## *** Check if we have the right format for plot.data
  if (class(plot.data) == 'data.frame'){
    if(sum(grepl('Name', names(plot.data))) != 1){
      stop("plot.data should contain 'Name' for complexes ... ")
    } else{
      if(is.null(cutoff.all)){
        print('Using precision cutoffs directly from file')
      }
      else if(dim(plot.data)[2] != (length(cutoff.all) + 1)){
        stop("plot.data and cutoff.all doesn't match ...")
      }
    }
  } else{
    stop("A data.frame with 'Name' and data at multiple precisions are expected ...")
  }
  
  if (!is.null(list.of.complexes.to.show)){
    if (num.complex.to.show != length(list.of.complexes.to.show) ){
      stop("Size of number of complexes to show and list of complexes doesn't match ...")
    }
  }
  
  # Remove duplicated data if not already done
  plot.data <- plot.data[!duplicated(plot.data$Name), ]
  
  ## *** Separate the Names and the TP numbers per precision
  cont_stepwise_anno <- plot.data$Name # Save the complex names (1824 * 1)
  cont_stepwise_mat <- plot.data[,-1] # 1824 * 38
  
  ## *** Remove precisions with less than min.pairs (10) pairs 
  tmp_TP <- apply(cont_stepwise_mat, 2, sum) # Summing the number of pairs at each cutoff
  Precision_ind <- (tmp_TP >= min.pairs) # Using 10 TP as the default parameter
  cont_stepwise_mat <- cont_stepwise_mat[,Precision_ind]
  
  # TODO: automatically get this from the data (plot.data)
  if (is.null(cutoff.all)){
    tmp <- names(cont_stepwise_mat)
    y <- as.numeric(substr(tmp, 11, max(sapply(tmp, nchar))))
  } else {y <- cutoff.all[Precision_ind] }
  
  # Find the fraction of TP contributed per complex (at a certain precision)
  # https://haky-functions.blogspot.com/2006/11/repmat-function-matlab.html
  x <- apply(cont_stepwise_mat, 2, sum)
  mx = dim(cont_stepwise_mat)[1] # 1824
  nx = dim(cont_stepwise_mat)[2] # 33
  tmp <- matrix(x, nrow = mx, ncol = nx, byrow = T)
  relx <- cont_stepwise_mat / tmp # 1824 * 33
  
  # Alternate version: Take out the background precision (not doing here)
  # x <- relx[, -1] # 1824 * 32 
  # y <- y[-1]
  x <- relx # 1824 * 32
  row.names(x) <- cont_stepwise_anno
  
  
  ## *** Arrange the data from smallest to largest contribution and 
  # ** Ranking: Use mean contributions from all precisions >= min.precision.cutoff (default 0.5) to rank the complexes
  ind.for.mean <- which(y >= min.precision.cutoff) # x and cutoff.all
  
  if (length(ind.for.mean) < 3){
    # If we don't have precision over 0.5, or very few of them, use the mid precision
    ind.for.mean <- which(y >= ((max(y)-min(y)) / 2) )
  }
  
  tmp.x <- x[,-1] # Not considering the background
  ind.for.mean <- ind.for.mean - 1 # Adjusting the indices  
  
  # Old ways: 1 (Rank by uniform representation of contribution (approximately) - Not used)
  # tmp.x <- x[,-1] # Not considering the background
  # ind.for.mean <- c(1) # Minimum (0.1)
  # ind.for.mean <- c(17) # 0.5
  # spaced.ind <- round(seq(from = ind.for.mean + 1, to = dim(tmp.x)[2]-1, length.out = min(8,dim(tmp.x)[2]-2)))
  # ind.for.mean <- c(ind.for.mean, spaced.ind) # min(8, remaining)
  # ind.for.mean <- c(ind.for.mean, dim(tmp.x)[2]) # Maximum (max precision with TP >= 10)
  
  # Old ways: 2 (Sort by mean across top 10 precisions)
  # tmp.x <- x
  # ind.for.mean <- (dim(x)[2]-9):dim(x)[2]
  
  # ** Take the bottom (largest contributions) num.complex.to.show (default 10)
  if (is.null(list.of.complexes.to.show)){
    a <- order(apply(tmp.x[,ind.for.mean], 1, mean))
    lx <- a[(length(a) - (num.complex.to.show - 1) ) : length(a)] 
    x <- x[lx, ]
  } else{ # If a list of complex is provided, use that!
    x <- x[list.of.complexes.to.show, ]
    if (!is.null(alternative.names)){ # Replace names with alternatives (if provided)
      row.names(x) <- alternative.names
    }
    
    x <- x[seq(dim(x)[1],1),] # Need to revert
    ccol <- ccol[seq(length(ccol),1)]
    
    # If there is NA
    non.na.ind <- !is.na(x[,1])
    x <- x[non.na.ind,] # If there is NA due to complex name not found!
    ccol <- ccol[non.na.ind]
  }
  
  ## *** Settle colors for top 10 complexes. If 10 colors are not provided, 
  #  we use a red to blue scale
  #  https://www.datanovia.com/en/blog/top-r-color-palettes-to-know-for-great-data-visualization/
  #  
  if (is.null(ccol)){
    # ccol <- c(colorRampPalette(colors = c("#bb2003","#3e7acf"))(dim(x)[1]))
    ccol <- c(colorRampPalette(colors = c("red","blue"))(dim(x)[1]))
  } else{
    if(length(ccol) < dim(x)[1]){
      warning ('Number of complex to show and number of colors provide do not match!! Using default coloring ...')
      ccol <- c(colorRampPalette(colors = c("red","blue"))(dim(x)[1]))
    }
  }
  
  ## *** Add small (rest of the) complexes
  xb <- x
  xb <- rbind(1 - apply(xb, 2, sum), xb) # Contribution of other complexes not in xb (14 * 37 now)
  row.names(xb)[1] <- "others"
  ccol <- c("#d9d9d9", ccol) # Gray (for rest)
  
  ## *** Data for Muller Plot
  x <- xb#[,-1]
  x1 <- x; x2 <- x
  
  for(i in 1:dim(x)[1]) {
    if(i == 1) {
      x1[i,] <- rep(0, dim(x)[2])
      x2[i,] <- x[1,]
    } else if(i == 2) {
      x1[i,] <- x[1,]
    } else {
      x1[i,] <- apply(x[1:(i - 1),], 2, sum) # Cumulative but doesn't include the last element; will be equal to the penultimate column of x2
    }
    
    if(i > 1) {
      x2[i,] <- apply(x[1:i,], 2, sum) # Cumulatively add complex contribution at each precision (the last value will be 1)
    }
  }
  
  # Save the figure
  if (save.figure){
    if(outfile.type == 'png'){
      png(paste0(outfile.name, ".png"), width = 4, height = 4, units="in", res = 300)
    }else{
      pdf(file = paste0(outfile.name, ".pdf"), width = 4, height = 4, useDingbats = F)  
    }
  }
  
  if (show.legend == TRUE){ # If we have to print legends (complex names)
    # Create a layout
    nf <- layout(matrix(c(1,2), 2,1), c(5,5), c(3,2), TRUE) #layout.show(nf)
    par(mar = c(4,4,1,1)) # bottom, left, top, right (order of margin)
    plot(as.double(x1[1,]), y, type = "l", xlim = c(0,1), ylim = y.lim, col = "white", bty = "n", las = 1, xlab = fig.labs[1], ylab = fig.labs[2], main = fig.title, cex.lab = 1.4, cex.axis = 1.4)
    
    for(i in 1:dim(x1)[1]) {
      polygon(c(x1[i,], rev(x2[i,])), c(y, rev(y)), col = ccol[i], border = "white")
    }
    
    # Restricting each complex to 55 chars so that it doesn't overflow
    legend.complex <- unlist(lapply(rev(row.names(x1)), substr, 1, 55)) 
    
    par(mar = c(1,1,0,1))
    # Make a null plot for legend
    plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1) 
    # Plot the legend (reverse so the best complex comes to the top)
    legend("center", legend = legend.complex, fill = rev(ccol), cex = 0.6, bty='n')
    
  } else{ # Just show the contribution plot (without any names)
    plot(as.double(x1[1,]), y, type = "l", xlim = c(0,1), ylim = y.lim, col = "white", bty = "n", las = 1, xlab = fig.labs[1], ylab = fig.labs[2], main = fig.title, cex.lab = 1.4, cex.axis = 1.4)
    
    for(i in 1:dim(x1)[1]) {
      polygon(c(x1[i,], rev(x2[i,])), c(y, rev(y)), col = ccol[i], border = "white")
    }
  }
  
  if (save.figure){
    dev.off()
  }
  
  return (rev(row.names(x1)))
}



#' Plot Category level PR Curve
#'
#' @param data Stepwise Contribution (numeric) for each complexes
#' @param anno Annotation (names) of the complexes in data
#' @param percent_th Normalized (by number of pairs) contribution for a complex to be considered as covered
#' @param tp_th Minimum number of TPs for a complex to be considered as covered
#' @param excludeBG Remove the background Precision? T/F
#' @param out_TP If false, plots recall instead of TP
#'
cTP_func <- function(data, anno, complexes,
                     percent_th = 0.05, tp_th = 1, excludeBG = T, out_TP = T) {
  
  # Replace modified names (for duplicate names that were modified before!)
  # Shouldn't happen anymore as we are removing the duplicates? But should we go back and modify the same name complexes?
  x <- which(anno %in% names(complexes) == F) 
  for(i in x) {
    anno[i] <- substr(anno[[i]], 1, nchar(anno[[i]]) - 2)
  }
  
  x <- data # Stepwise contribution per precision
  for(i in 1:dim(data)[1]) {
    # print(data[i,])
    # print(choose(length(complexes[[anno[[i]]]]), 2))
    
    # Normalizing complex contribution by number of possible gene pairs (nC2)
    x[i,] <- data[i,] / choose(length(complexes[[anno[[i]]]]), 2) 
  }
  
  # Applying a TP cutoff and a percent of complex contribution cutoff
  x <- c(apply(data >= tp_th & x > percent_th, 2, sum), 1)
  
  if(out_TP == FALSE) {
    x <- x / max(x) #recall
  }
  
  if(excludeBG) {
    x <- x[-which(x == max(x))]
  }
  
  return (x)
}


#' Separate the contribution of signal into individual entries (a complex, pathway, GO Biological Process etc.)
#'
#' @param data_complex Original input for the co-annotation standard
#' @param pr.stepwise Stepwise Contribution (TP) per complex at different precisions
#' @param legend.names legends for different curves
#' @param ccol colors for each of the PR curves
#' @param thresholds a vector of two thresholds to filter complexes (1) No of TP, (2) Percent of TP pairs captured.
#' @param save.figure set to TRUE to save the figure (fig.title is used as the name)
#' @param outfile.name the name of the output file(figure)
#' @param outfile.type type of figure to save - 'pdf' (default) or 'png'
#' 
#' @examples
#' data('data_complex', package = 'FLEX')
#' pr_full <- read.table(paste0('Stepwise_contribution_Complex_', datasets[i] ,'.txt'), stringsAsFactors=FALSE, sep = "\t", header = T)
#' pr_removed <- read.table(paste0('Stepwise_contribution_Complex_Removal_', datasets[i] ,'.txt'), stringsAsFactors=FALSE, sep = "\t", header = T)
#' pr.stepwise <- list(full = list(data = pr_full))
#' pr.stepwise <- append(pr.stepwise, list(removed = list(data = pr_removed)))
#' PlotCategoryPR(data_complex, pr.stepwise, thresholds = c(1, 0.1), ccol = c('#252525', '#bd0026'), fig.labs = c('TP Complexes', 'Precision'))
#' 
#' @export
#' 
PlotCategoryPR <- function(data_complex, pr.stepwise, thresholds = NULL, 
                           legend.names = NULL, legend.ltype = NULL, 
                           ccol = NULL, 
                           fig.labs = c('Category TP', 'Precision'),
                           save.figure = FALSE, 
                           outfile.name = 'test_category_PR', outfile.type = 'pdf'){
  
  # Remove complexes (say top complex, top 3, top 5, top 10 etc.) and generate the curve shown on the picture
  coreComplex_members <- strsplit(data_complex$Genes, "[;]") # ';' would work just fine
  names(coreComplex_members) <- data_complex$Name
  
  # Save figure
  if (save.figure){
    if(outfile.type == 'png'){
      png(paste0(outfile.name, ".png"), width = 2.5, height = 3, units="in", res = 300)
    }else{
      pdf(file = paste0(outfile.name, ".pdf"))
    }
  }
  
  # Get the xlim for TP plot (maximum among all standard)
  x.lim <- 1
  
  # First get the maximum TP among the standards to fix the x.lim
  for (i in 1 : length(pr.stepwise)) {
    pr_contri = pr.stepwise[[i]]$data
    pr_contri_anno <- pr_contri$Name
    pr_contri <- pr_contri[,-1] # Remove complex Name, and keep the data
    
    ## Old way
    # y <- pr.stepwise[[i]]$cutoffs
    
    ## New way (calculation of y) - now we don't need to send cutoffs
    tmp <- names(pr_contri)
    y <- as.numeric(substr(tmp, 11, max(sapply(tmp, nchar))))
    y <- y[-1] # no background is shown
    y <- c(y,1) # A 1 is added! Why!
    
    if (is.null(thresholds)){
      x <- cTP_func(data = pr_contri, anno = pr_contri_anno,
                    complexes = coreComplex_members, percent_th = .3)
    } else{
      x <- cTP_func(data = pr_contri, anno = pr_contri_anno,
                    complexes = coreComplex_members, 
                    tp_th = thresholds[1], percent_th = thresholds[2])
    } # else
    
    print(max(x))
    
    if (max(x) > x.lim){
      x.lim <- max(x)
    }
    
    # Save the data so that we don't have to run twice!
    if (i == 1){
      data_to_out <- list(data = x)
    } else{
      data_to_out <- append(data_to_out, list(data = x))
    }
    
  } 
  x.lim

  # Now do the plotting
  for (i in 1 : length(pr.stepwise)) {
    
    x <- data_to_out[[i]]
    
    if (i == 1){
      plot(x, y, xlab = fig.labs[1], ylab = fig.labs[2], bty = "n", las=1, 
           ylim = c(0,max(y)), xlim = c(1,x.lim),
           pch = 16, 
           type = "l", log = "x", lwd = 2,
           cex.lab = 1.4, cex.main = 1.4, cex.axis = 1.4,
           col = ccol[i], lty = legend.ltype[i])
      
    }else{
      lines(x, y, col = ccol[i], lwd = 2,)
    }
  }
  
  if(!is.null(legend.names)){
    legend("topright", legend = legend.names, fill = ccol, 
           cex = 1, bty = "n", text.col = "black", horiz = F)
  }

if (save.figure){
  dev.off()
}

}
