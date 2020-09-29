## ============================ ##
## Test code for FLEX R Package
## ============================ ##

## Load necessary fiels

# Way 1: Source the necessary files directly
# only works when we source this file
# test.dir <- dirname(sys.frame(1)$ofile)
# src.dir <- gsub('tests', 'R', test.dir)
# setwd(src.dir)
# setwd(test.dir)

# This directories are subject to change
setwd('/project/chadm/Mahfuz/CRISPR/2_HAP1/GIN_Analysis/FLEX/R/FLEX_R/R')
source('Utility.R')
source('GoldStandards.R')
source('AssembleDataForPR.R')
source('PerfCurve.R')
source('PostProcessing.R')
source('Plots.R')

setwd('/project/chadm/Mahfuz/CRISPR/2_HAP1/GIN_Analysis/FLEX/R/FLEX/tests')

# Way 2 (use devtools to source R files from FLEX)
require('devtools')
# load_all('/home/mahfuz/Desktop/CRISPR/FLEX/R/FLEX_R')
load_all('/project/chadm/Mahfuz/CRISPR/2_HAP1/GIN_Analysis/FLEX/R/FLEX_R/R')

# Way 3: load the FLEX library
# install('FLEX') # From the directory where it resides
# detach('package:FLEX', unload = TRUE)
# library(FLEX)


## ================================================ ##
#  Demo of the FLEX package using Profile Similarity
## ================================================ ##

# ---------------------------------------------
# 1. Read in the data included with the package
# ---------------------------------------------
data('data_complex', package = 'FLEX')

# or load it directly
# setwd('/project/chadm/Mahfuz/CRISPR/2_HAP1/GIN_Analysis/FLEX/R/FLEX/data/')
# load('data_complex.rda')

# ---------------------------------------------
# 2. Create (or read) the co-annotation data
# ---------------------------------------------
file_name <- '/project/chadm/Mahfuz/CRISPR/2_HAP1/GIN_Analysis/FLEX/R/CORUM_CA.Rdata'
file_name <- '/home/mahfuz/Desktop/CRISPR/FLEX/R/Tests/Package_Test/CORUM_CA.Rdata'
data.ca <- MakeCoAnnotationFromGeneSymbols(data_standard = data_complex, 
                                           overlap_length = 1, 
                                           file_location = file_name)

# ---------------------------------------------
# 3. Read in the interaction/dependency data
# ---------------------------------------------

# We are going to use the sample data included with the package here
data('data_depmap_19Q2_test', package = 'FLEX')

# Alternative (load data from text file)
# file.int <- 'Ceres_score_19Q2_with_depmap_id.txt' # tab delimited 
# data.interaction <- GetInteractionData(file.int)

# ---------------------------------------------
# 4. Associate the pairwise scores to co-annotation
# ---------------------------------------------
if (!file.exists('Complex_Depmap_19Q2_small.Rdata')){
  
  # Way 1: Using the dependency data directly
  Complex.DepMap.19Q2.set1 <- CalculatePredictionAndTrueOnLibraryProfiles (data.ca, data.interaction.set1)
  
  # Way 2: Using pre-calculated pairwise correlation matrix
  # start_time <- Sys.time()
  # pairwise.correlation <- cor(t(data.interaction), use = 'pairwise.complete.obs', method = 'pearson')
  # end_time <- Sys.time()
  # end_time - start_time
  # Complex.DepMap.19Q2 <- CalculatePredictionAndTrueOnLibraryProfiles (data.ca, pairwise.correlation)
  
  save(Complex.DepMap.19Q2.set1, file = 'Complex.DepMap.19Q2.set1.Rdata')
}else{
  load('Complex.DepMap.19Q2.set1.Rdata')
}

# Same for another dataset (TO make two PR Curves)
Complex.DepMap.19Q2.set2 <- CalculatePredictionAndTrueOnLibraryProfiles (data.ca, data.interaction.set2)
save(Complex.DepMap.19Q2.set2, file = 'Complex.DepMap.19Q2.set2.Rdata')

# For cont str plot test (with smaller precisions)
# data.interaction.pca <- GetInteractionData('normalized_pca.tsv')
# Complex.DepMap.19Q2.pca <- CalculatePredictionAndTrueOnLibraryProfiles (data.ca, data.interaction.pca)
# save(Complex.DepMap.19Q2.pca, file = 'Complex.DepMap.19Q2.pca.Rdata')

# ---------------------------------------------
# 5. Plot global PR Curves
# ---------------------------------------------
load('Complex.DepMap.19Q2.set1.Rdata')
load('Complex.DepMap.19Q2.set2.Rdata')

pred.ca <- list(out_19Q2_1 = list(true = Complex.DepMap.19Q2.set1$true, 
                                predicted = Complex.DepMap.19Q2.set1$predicted))
pred.ca <- append(pred.ca, list(out_19Q2_2 = list(true = Complex.DepMap.19Q2.set2$true, 
                                                predicted = Complex.DepMap.19Q2.set2$predicted)))
# pred.ca <- append(pred.ca, list(out_19Q2_pca = list(true = Complex.DepMap.19Q2.pca$true, predicted = Complex.DepMap.19Q2.pca$predicted)))

PlotPRSimilarity (pred.ca, subsample = TRUE, type.plot = 'log',
                  fig.title = 'DepMap 19Q2 subsets', 
                  fig.labs = c('TP', 'Precision'), legend.names = c('Set 1', 'Set 2'), 
                  legend.color = c('#de2d26', '#3182bd'), save.figure = FALSE)

# PlotPRSimilarity (pred.ca, fig.title = 'DepMap 19Q2 subsets',
#                   fig.labs = c('TP', 'Precision'), legend.names = c('Set 1', 'Set 2'),
#                   legend.color = c('#de2d26', '#3182bd'), save.figure = TRUE,
#                   outfile.name = 'DepMap_19Q2', outfile.type = 'pdf')

# PlotPRSimilarity (pred.ca, subsample = TRUE, type.plot = 'log', legend.names = c('Set 1', 'Set 2', 'PCA'), legend.color = c('#de2d26', '#3182bd'), save.figure = FALSE)

# ---------------------------------------------
# 6. Individual AUPRC (Contribution Scatter)
# ---------------------------------------------

entity.matrix <- as.data.frame(Complex.DepMap.19Q2.set1, stringsAsFactors = FALSE)
data.AUPRC <- GetAreaUnderPRCurveForEntities (summary.standard = data_complex, data.standard = data.ca, entity.matrix = entity.matrix)
#write.table(data.AUPRC, 'Complex_AUPRC_DepMap_19Q2_set1.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#data.AUPRC <- read.table('Complex_AUPRC_DepMap_19Q2_set1.txt', stringsAsFactors=FALSE, sep = "\t", header = T, quote = '')
plot.data <- data.AUPRC[!duplicated(data.AUPRC$Name), ]

PlotContributionScatter (plot.data, fig.labs = c('AUPRC', 'Complex size'), show.text = FALSE, save.figure = FALSE)

# Don't advise using show.text = TRUE, looks like a mess
#PlotContributionScatter (plot.data, fig.labs = c('AUPRC', 'Complex size'), show.text = TRUE, save.figure = TRUE, outfile.type = 'pdf', outfile.name = 'test_scatter')

# ---------------------------------------------
# 7. Contribution Structure Plot
# ---------------------------------------------
Pairs.in.data <- data.frame(true = Complex.DepMap.19Q2.set1$true, 
                            predicted = Complex.DepMap.19Q2.set1$predicted, 
                            ID = Complex.DepMap.19Q2.set1$ID, stringsAsFactors = FALSE)

#Pairs.in.data <- data.frame(true = Complex.DepMap.19Q2.pca$true, predicted = Complex.DepMap.19Q2.pca$predicted, ID = Complex.DepMap.19Q2.pca$ID, stringsAsFactors = FALSE)

# Find the Precision cutoffs to use (be careful about this)
out_19Q2 <- GenerateDataForPerfCurve(value.predicted = Complex.DepMap.19Q2.set1$predicted, 
                                        value.true = Complex.DepMap.19Q2.set1$true, 
                                        x.axis = 'TP', y.axis = 'precision')

precision_cutoffs <- c(out_19Q2$y[length(out_19Q2$y)], seq(0.1, max(out_19Q2$y), 0.025)) # [bgd_precision 0.1, 0.125, 0.15, ..., max_precision]
# precision_cutoffs <- c(out_19Q2$y[length(out_19Q2$y)], seq(0.1, max(out_19Q2$y), length = 20)) # [bgd_precision 0.1 18 in between 1.0]
precision_cutoffs <- round(precision_cutoffs, 3)

# Generate the stepwise contribution of complexes
output.stepwise.contribution <- GetStepwiseContributionOfEntities(Pairs.in.data = Pairs.in.data, cutoff.all = precision_cutoffs, summary.standard = data_complex)
#write.table(output.stepwise.contribution, 'Contribution_of_complexes_stepwise_DepMap_19Q2.txt', sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

#output.stepwise.contribution <- read.table('Contribution_of_complexes_stepwise_DepMap_19Q2.txt', stringsAsFactors=FALSE, sep = "\t", header = T)
output.stepwise.contribution <- output.stepwise.contribution[!duplicated(output.stepwise.contribution$Name), ]


PlotContributionStructure(plot.data = output.stepwise.contribution, cutoff.all = precision_cutoffs, min.pairs = 10, fig.title = 'Contriubtion Structure Plot DepMap 19Q2', save.figure = FALSE)

# PlotContributionStructure(plot.data = output.stepwise.contribution, cutoff.all = precision_cutoffs, min.pairs = 10, fig.title = 'Contriubtion Structure Plot DepMap 19Q2', save.figure = TRUE, outfile.name = 'Str_Depmap_19Q2', outfile.type = 'pdf')

# Same but using custom colors for 10 complexes
#PlotContributionStructure(plot.data = output.stepwise.contribution, cutoff.all = precision_cutoffs, min.pairs = 10, save.figure = FALSE, ccol = c('#084594', '#4292c6', '#c6dbef', '#005a32', '#41ab5d', '#a1d99b', '#99000d', '#ef3b2c', '#fc9272', '#4a1486'), fig.labs = c('Fraction TP', 'Precision'), outfile.name = 'Contriubtion Structure Plot DepMap 19Q2')

#PlotContributionStructure(plot.data = output.stepwise.contribution, cutoff.all = precision_cutoffs, min.pairs = 10, save.figure = FALSE, ccol = rainbow(10)) # rainbow, heat.colors, terrain.colors, topo.colors, cm.colors

# ---------------------------------------------
# 8. Complex Removal analysis
# ---------------------------------------------

# Remove top 10 complexes (according to AUPRC values)
data.AUPRC <- read.table('Complex_AUPRC_DepMap_19Q2_set1.txt', stringsAsFactors=FALSE, sep = "\t", header = T, quote = '')
data.AUPRC.nonDuplicated <- data.AUPRC[!duplicated(data.AUPRC$Name), ]

# sort by AUPRC (should already be though)
data.AUPRC.nonDuplicated <- data.AUPRC[order(data.AUPRC.nonDuplicated$AUPRC, decreasing = F),]

# top.complex.ids <- as.character(data.AUPRC.nonDuplicated$ID[1:10])
entity.matrix <- as.data.frame(Complex.DepMap.19Q2.set1, stringsAsFactors = FALSE)
top.complex.ids <- c('320')
entity.matrix.top10.removed <- getSubsetOfCoAnnRemoveIDs(data.standard = entity.matrix, ids = top.complex.ids, replace = FALSE)


# PR Curve comparison
pred.ca <- list(out_orig = list(true = Complex.DepMap.19Q2.set1$true, predicted = Complex.DepMap.19Q2.set1$predicted))
pred.ca <- append(pred.ca, list(out_cmplx_rem_10 = list(true = entity.matrix.top10.removed$true, predicted = entity.matrix.top10.removed$predicted)))

PlotPRSimilarity (pred.ca, fig.title = 'Complex Removal Comparison', fig.labs = c('TP', 'Precision'), legend.names = c('Original', '55S_removed'), legend.color = c('#de2d26', '#3182bd'), save.figure = FALSE)


# ---------------------------------------------
#  9. Category Plot
# ---------------------------------------------
setwd('/home/mahfuz/Desktop/CRISPR/FLEX/R/Tests/Package_Test/Vignette_test')

## Stepwise contribution (Data inside Vignette directory)
pr_contri <- read.table('Contribution_of_complexes_stepwise_19Q2.txt', stringsAsFactors=FALSE, sep = "\t", header = T)
pr_contri_noETC <- read.table("Contribution_of_complexes_stepwise_19Q2_ETC1_mtRibo_ETCV_removal.txt", stringsAsFactors=FALSE, sep = "\t", header = T)
pr_contri_noAUChi <- read.table("Contribution_of_complexes_stepwise_19Q2_low_size_high_AUC_removal.txt", stringsAsFactors=FALSE, sep = "\t", header = T)

# The second parameter should match up with the precision values used to generate pr_contri
pr.stepwise <- list(first = list(data = pr_contri, cutoffs = c(seq(.1,1,0.025), 1)))
pr.stepwise <- append(pr.stepwise, list(second = list(data = pr_contri_noETC, cutoffs = c(seq(.1,1,0.025), 1))))
pr.stepwise <- append(pr.stepwise, list(third = list(data = pr_contri_noAUChi, cutoffs = c(seq(.1,1,0.025), 1))))

PlotCategoryPR(data_complex, pr.stepwise, thresholds = c(1, 0.3), ccol = c('#252525', '#fb6a4a', '#74c476'), legend.names = c('Full', 'ETC removed', 'high AUC, low size removed'))



## Following analysis are for testing the package
if (FALSE){
  ## =========================== ##
  #  Direct Interaction 
  ## =========================== ##
  
  # filename.co.ann <- '/project/chadm/Mahfuz/Database/Custom/Pathway/Pathway_CO_Annotation_MsigDB_Symbol_Pathway_ID.txt' # Pathway
  filename.co.ann <- '/project/chadm/Mahfuz/Database/Custom/Complex/v3.0/CORUM_Human_CO_Annotation_Symbol_ID.txt';
  data.standard <- GetCoAnnotationData(filename.co.ann)
  
  # Interaction data: 2
  filename.interaction <- '/project/chadm/Mahfuz/CRISPR/2_HAP1/GIN_Analysis/Scoring_Manuscript/Data/qGI_20200212.txt'
  data.interaction <- GetInteractionData(filename.interaction, '', '')
  
  direct.pairing.qGI_20200212 <- CalculatePredictionAndTrueOnDirectInteraction (data.standard, data.interaction)
  save(direct.pairing.qGI_20200212, file = "DI_AUC_direct.pairing.qGI_20200212.RData")
  
  # Load the data for plots
  load ('/project/chadm/Mahfuz/CRISPR/2_HAP1/GIN_Analysis/BROAD_vs_TKOV3/For_Manuscript/Package/R/DI_AUC_direct.pairing.20200121.RData')
  
  # Plot co-annotation 
  PlotPRDirect (direct.pairing.qGI_20200212)
  
  
  ## Use PR similarity plot for direct interactions
  filename.co.ann <- '/project/chadm/Mahfuz/Database/Custom/Complex/v3.0/CORUM_Human_CO_Annotation_Symbol_ID.txt';
  data.standard <- GetCoAnnotationData(filename.co.ann)
  
  # First data
  filename.interaction <- '/project/chadm/Mahfuz/CRISPR/2_HAP1/GIN_Analysis/Scoring_Manuscript/Data/qGI_20200212.txt'
  filename.queries <- '/project/chadm/Mahfuz/CRISPR/2_HAP1/GIN_Analysis/Scoring_Manuscript/final_unique_screens_with_shifts_ID.txt'
  data.interaction <- GetInteractionData(filename.interaction, filename.queries = filename.queries)
  direct.pairing.qGI_20200212_small <- CalculatePredictionAndTrueOnDirectInteraction (data.standard, data.interaction)
  save(direct.pairing.qGI_20200212_small, file = "DI_AUC_direct.pairing.qGI_20200212_small_queries.RData")
  
  # Second data
  filename.interaction <- '/project/chadm/Mahfuz/CRISPR/2_HAP1/GIN_Analysis/Scoring_Manuscript/Data/qGI_noCorr_20200212.txt'
  filename.queries <- '/project/chadm/Mahfuz/CRISPR/2_HAP1/GIN_Analysis/Scoring_Manuscript/final_unique_screens_with_shifts_ID.txt'
  data.interaction <- GetInteractionData(filename.interaction, filename.queries = filename.queries)
  direct.pairing.qGI_20200212_noCorr_small <- CalculatePredictionAndTrueOnDirectInteraction (data.standard, data.interaction)
  save(direct.pairing.qGI_20200212_noCorr_small, file = "DI_AUC_direct.pairing.qGI_20200212_noCorr_small_queries.RData")
  
  
  # Make the data
  pred.ca <- list(small = list(true = direct.pairing.qGI_20200212_small$data$True, predicted = direct.pairing.qGI_20200212_small$data$Score))
  pred.ca <- append(pred.ca, list(noCorr_small = list(true = direct.pairing.qGI_20200212_noCorr_small$data$True, predicted = direct.pairing.qGI_20200212_noCorr_small$data$Score)))
  
  
  # Compare between two datasets (Neg-neg, or pos-pos)
  # To get good colors: http://colorbrewer2.org/#type=sequential&scheme=BuGn&n=6
  legend.names <- c('qGI Score', 'qGI Score noCorr')
  rgbConversion <- (function(x,y,z) rgb(x/255, y/255, z/255))
  fig.labs = c('TP', 'Precision')
  type.plot = 'log'
  
  # Positive
  fig.title <- 'Complex (Positive) Interaction'
  legend.col.pos <- c(rgbConversion(255, 222, 23), rgbConversion(255, 255, 178))
  PlotPRSimilarity (pred.ca = pred.ca, neg.to.pos = FALSE, fig.title = fig.title, legend.names = legend.names, legend.color = legend.col.pos, save.figure = FALSE, is.bgdline = TRUE)
  
  # Negative
  fig.title <- 'Complex (Negative) Interaction'
  legend.col.neg <- c(rgbConversion(0, 118, 239), rgbConversion(158, 202, 225))
  
  PlotPRSimilarity (pred.ca = pred.ca, neg.to.pos = TRUE, fig.title = fig.title, legend.names = legend.names, legend.color = legend.col.neg, save.figure = FALSE, is.bgdline = TRUE)
}


if (FALSE){
  ## PerfCurve Test (Compare with Matlab for GLS)
  load('/home/mahfuz/Desktop/CRISPR/FLEX/Analysis/Supplementary_Figures/S5/Contr_Str_Plot/PR_Similarity_Complex_GLS_weighted.RData')
  test <- GenerateDataForPerfCurve(value.predicted = Complex.GLS.weighted$predicted, 
                                   value.true = Complex.GLS.weighted$true, 
                                   x.axis = 'TP', y.axis = 'precision')
  
  pred.ca <- list(out_orig = list(true = Complex.GLS.weighted$true, predicted = Complex.GLS.weighted$predicted))
  
  PlotPRSimilarity (pred.ca, subsample = TRUE, fig.title = 'Complex Removal Comparison', fig.labs = c('TP', 'Precision'), legend.names = c('Original', '55S_removed'), legend.color = c('#de2d26', '#3182bd'), save.figure = FALSE)
}

