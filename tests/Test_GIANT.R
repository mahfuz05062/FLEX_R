## ============================ ##
## Test code for FLEX R Package
## ============================ ##

require(devtools)
# load_all('/home/mahfuz/Desktop/CRISPR/FLEX/R/FLEX/') # Roughly simulates library('FLEX')
load_all('/project/chadm/Mahfuz/CRISPR/2_HAP1/GIN_Analysis/FLEX/R/FLEX_R/')

# scp rahma118@mcclintock.cs.umn.edu:/project/chadm/Mahfuz/CRISPR/2_HAP1/GIN_Analysis/FLEX/R/FLEX/R/*.R /home/mahfuz/Desktop/CRISPR/FLEX/R/FLEX/R/
# scp /home/mahfuz/Desktop/CRISPR/FLEX/R/FLEX/R/*.R rahma118@mcclintock.cs.umn.edu:/project/chadm/Mahfuz/CRISPR/2_HAP1/GIN_Analysis/FLEX/R/FLEX/R/ 
# scp /home/mahfuz/Desktop/CRISPR/FLEX/R/FLEX/tests/*.R rahma118@mcclintock.cs.umn.edu:/project/chadm/Mahfuz/CRISPR/2_HAP1/GIN_Analysis/FLEX/R/FLEX/tests/ 


# If for example, devtools is not installed
if (FALSE){
  setwd('/project/chadm/Mahfuz/CRISPR/2_HAP1/GIN_Analysis/FLEX/R/FLEX/R')
  source('Utility.R')
  source('GoldStandards.R')
  source('AsembleDataForPR.R')
  source('PerfCurve.R')
  source('PostProcessing.R')
  source('Plots.R')
  setwd('/project/chadm/Mahfuz/CRISPR/2_HAP1/GIN_Analysis/FLEX/R/FLEX/data')
}

# 1. Read in the data included with the package
# ---------------------------------------------
data('data_pathway', package = 'FLEX')


# 2. Create (or read once created) the co-annotation data
# -------------------------------------------------------
# file_name <- '/home/mahfuz/Desktop/CRISPR/FLEX/R/Tests/Package_Test/Pathway_CA.Rdata'
file_name <- '/project/chadm/Mahfuz/CRISPR/2_HAP1/GIN_Analysis/FLEX/R/Pathway_Test/Pathway_CA.Rdata'
data.ca <- MakeCoAnnotationFromGeneSymbols(data_standard = data_pathway, 
                                           overlap_length = 1, 
                                           file_location = file_name)

# 3. Read in the interaction/dependency data
# ---------------------------------------------
# We are going to use the sample data included with the package here
data('data_depmap_19Q2_test', package = 'FLEX')

# Alternative (load data from a text file)
# file.int <- 'Ceres_score_19Q2_with_depmap_id.txt' # tab delimited 
# data.interaction <- GetInteractionData(file.int)


# 4. Associate the pairwise scores to co-annotation
# ---------------------------------------------
setwd('/project/chadm/Mahfuz/CRISPR/2_HAP1/GIN_Analysis/FLEX/R/Pathway_Test')

Pathway.DepMap.19Q2.set1 <- CalculatePredictionAndTrueOnLibraryProfiles (data.ca, data.interaction.set1)
save(Pathway.DepMap.19Q2.set1, file = 'Pathway.DepMap.19Q2.set1.Rdata')

Pathway.DepMap.19Q2.set2 <- CalculatePredictionAndTrueOnLibraryProfiles (data.ca, data.interaction.set2)
save(Pathway.DepMap.19Q2.set2, file = 'Pathway.DepMap.19Q2.set2.Rdata')


# 5. Plot global PR Curves
# ---------------------------------------------
load('Pathway.DepMap.19Q2.set1.Rdata')
load('Pathway.DepMap.19Q2.set2.Rdata')

pred.ca <- list(out_19Q2_1 = list(true = Pathway.DepMap.19Q2.set1$true, 
                                  predicted = Pathway.DepMap.19Q2.set1$predicted))
pred.ca <- append(pred.ca, list(out_19Q2_2 = list(true = Pathway.DepMap.19Q2.set2$true, 
                                                  predicted = Pathway.DepMap.19Q2.set2$predicted)))

PlotPRSimilarity (pred.ca, subsample = TRUE, type.plot = 'log',
                  fig.title = 'DepMap 19Q2 subsets', 
                  fig.labs = c('TP', 'Precision'), legend.names = c('Set 1', 'Set 2'), 
                  legend.color = c('#de2d26', '#3182bd'), save.figure = TRUE,
                  outfile.name = 'DepMap_19Q2_Pathway', outfile.type = 'pdf')

# 6. Individual AUPRC (Contribution Scatter)
# ---------------------------------------------
entity.matrix <- as.data.frame(Pathway.DepMap.19Q2.set1, stringsAsFactors = FALSE)
data.AUPRC <- GetAreaUnderPRCurveForEntities (summary.standard = data_pathway, data.standard = data.ca, entity.matrix = entity.matrix)
write.table(data.AUPRC, 'Pathway_AUPRC_DepMap_19Q2_set1.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#data.AUPRC <- read.table('Pathway_AUPRC_DepMap_19Q2_set1.txt', stringsAsFactors=FALSE, sep = "\t", header = T, quote = '')
plot.data <- data.AUPRC[!duplicated(data.AUPRC$Name), ]

# Don't advise using show.text = TRUE, looks like a mess
PlotContributionScatter (plot.data, length.cutoff = 10, AUPRC.cutoff = 0.1, fig.labs = c('AUPRC', 'Pathway size'), show.text = FALSE, save.figure = TRUE, outfile.name = 'Pathway_Scatter_Plot')

# scp rahma118@mcclintock.cs.umn.edu:/project/chadm/Mahfuz/CRISPR/2_HAP1/GIN_Analysis/FLEX/R/Pathway_Test/*.txt /home/mahfuz/Desktop/CRISPR/FLEX/R/Tests/Package_Test



## This is super slow for pathway, so commented out for now!
if (FALSE){
  # 7. Contribution Structure Plot
  # ---------------------------------------------
  Pairs.in.data <- data.frame(true = Pathway.DepMap.19Q2.set1$true, 
                              predicted = Pathway.DepMap.19Q2.set1$predicted, 
                              ID = Pathway.DepMap.19Q2.set1$ID, stringsAsFactors = FALSE)
  
  # Find the Precision cutoffs to use (be careful about this)
  out_19Q2 <- GenerateDataForPerfCurve(value.predicted = Pathway.DepMap.19Q2.set1$predicted, 
                                       value.true = Pathway.DepMap.19Q2.set1$true, 
                                       x.axis = 'TP', y.axis = 'precision')
  
  # precision_cutoffs <- c(out_19Q2$y[length(out_19Q2$y)], seq(0.1, max(out_19Q2$y), 0.025)) 
  precision_cutoffs <- c(out_19Q2$y[length(out_19Q2$y)], seq(0.1, max(out_19Q2$y), length = 10)) # [bgd_precision 0.1 18 in between 1.0]
  precision_cutoffs <- round(precision_cutoffs, 3)
  
  # Generate the stepwise contribution of complexes
  output.stepwise.contribution <- GetStepwiseContributionOfEntities(Pairs.in.data = Pairs.in.data, cutoff.all = precision_cutoffs, summary.standard = data_pathway)
  write.table(output.stepwise.contribution, 'Contribution_of_pathways_stepwise_DepMap_19Q2.txt', sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  #output.stepwise.contribution <- read.table('Contribution_of_complexes_stepwise_DepMap_19Q2.txt', stringsAsFactors=FALSE, sep = "\t", header = T)
  output.stepwise.contribution <- output.stepwise.contribution[!duplicated(output.stepwise.contribution$Name), ]
  
  
  PlotContributionStructure(plot.data = output.stepwise.contribution, cutoff.all = precision_cutoffs, min.pairs = 10, fig.title = 'Contriubtion Structure Plot DepMap 19Q2', save.figure = TRUE, outfile.name = 'Str_Pathway_Depmap_19Q2', outfile.type = 'pdf')
}
