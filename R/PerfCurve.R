## Performance Curve with Co Annotation: to generate Co-Annotation PR Curves
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


#' Generate the values for performance curve (Precision and Recall for example) to plot on x and y axis
#'
#' @param value.predicted Predicted value / Score (profile similarity / direct interaction)
#' @param value.true Actual value (1/0 co-annotation for example)
#' @param neg.to.pos assign TRUE to sort from neg to positive (pos -> neg is the default)
#' @param x-axis what is the x-axis value you want? (TP, FP, TN, FN, precision, FPR, recall/sensitivity)
#' @param y-axis what is the y-axis value you want? (TP, FP, TN, FN, precision, FPR, recall/sensitivity)
#'
#' @return PR -> (x: recall, y: precision), ROC -> (x: FPR or 1 - specificity, and y: recall/TPR)
#' @examples
#' PerfCurve <- function(value.predicted, value.true, neg.to.pos = FALSE, type = 'PR', x.axis = 'sensitivity', y.axis = 'precision')
#' @export
#' 
GenerateDataForPerfCurve <- function(value.predicted, value.true, neg.to.pos = FALSE, x.axis = 'sensitivity', y.axis = 'precision'){
  
  # Sorting direction
  if (neg.to.pos == FALSE){
    indices <- order(value.predicted, decreasing = TRUE)
  } else{
    indices <- order(value.predicted)
  }
  value.true <- value.true[indices]
  value.predicted <- value.predicted[indices]
  
  # Calculate basic elements
  num.real.true <- sum(value.true)
  num.predicted.true <- 1 : length(value.true) # Predicted Positive <= Threshold
  
  TP <- cumsum(value.true)
  FP <- num.predicted.true - TP
  FN <- num.real.true - TP
  TN <- length(value.true) - (TP + FP + FN)
  
  # *** Get the indices of unique predicted values (last occurrence)
  # For cases when we have the same prediction for many elements
  unique.index <- which(!duplicated(value.predicted)) # Gives index of first occurrence
  unique.index <- c(unique.index - 1, length(TP)) # Now, it's the last occurrence.
  unique.index <- unique.index[-1] # The first element is 0 anyway!
  
  TP <- TP[unique.index]
  FP <- FP[unique.index]
  FN <- FN[unique.index]
  TN <- TN[unique.index]
  
  precision <- TP / (TP + FP) # PPV
  sensitivity <- TP / (TP + FN)  # Recall / TPR
  FPR <- FP / (FP + TN) # FPR 
  
  ## *** TODO: What about the cases when we have very few unique precision and/or sensitivity values?
  
  # Find which to return for x and y
  switch(x.axis, TP = {x = TP}, FP = {x = FP}, TN = {x = TN}, FN = {x = FN},
         precision = {x = precision}, 
         FPR = {x = FPR}, 
         sensitivity = {x = sensitivity}, recall = {x = sensitivity}, TPR = {x = sensitivity})
  switch(y.axis, TP = {y = TP}, FP = {y = FP}, TN = {y = TN}, FN = {y = FN},
         precision = {y = precision}, 
         FPR = {y = FPR}, 
         sensitivity = {y = sensitivity}, recall = {y = sensitivity}, TPR = {y = sensitivity})
  
  # area under curve: according to trapizoidal approximation (make sense for roc or pr curve only: unit area)
  # https://www.r-bloggers.com/calculating-auc-the-area-under-a-roc-curve/
  # area of trapezoid = 0.5 * h * (b1 + b2) # diff in x = h (height) # b1, b2 are values on the y axis surrounding a trapezoid
  auc = abs(0.5 * sum((x[-1] - x[-length(x)]) * (y[-1] + y[-length(y)]))) # same formula as used in matlab (and gives the same value)
  
  # Or this works too
  # dx <- c(diff(x), 0)
  # dy <- c(diff(y), 0)
  # sum(y * dx) + sum(dx * dy)/2
  
  return(list(x = x, y = y, auc = auc))
}