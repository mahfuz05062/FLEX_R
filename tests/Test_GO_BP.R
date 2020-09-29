## ============================ ##
## Test code for FLEX R Package
## ============================ ##

require(devtools)
# load_all('/home/mahfuz/Desktop/CRISPR/FLEX/R/FLEX/') # Roughly simulates library('FLEX')
load_all('/project/chadm/Mahfuz/CRISPR/2_HAP1/GIN_Analysis/FLEX/R/FLEX/') 

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
data('data_GO_BP', package = 'FLEX')
# load('data_GO_BP.rda')

# 2. Create (or read once created) the co-annotation data
# -------------------------------------------------------
file_name <- '/project/chadm/Mahfuz/CRISPR/2_HAP1/GIN_Analysis/FLEX/R/Pathway_Test/GO_BP_CA.Rdata'
data.ca <- MakeCoAnnotationFromGeneSymbols(data_standard = data_GO_BP, 
                                           overlap_length = 1, 
                                           file_location = file_name)


# 3. Read in the interaction/dependency data
# ---------------------------------------------
# We are going to use the sample data included with the package here
data('data_depmap_19Q2_test', package = 'FLEX')


# 4. Associate the pairwise scores to co-annotation
# ---------------------------------------------
setwd('/project/chadm/Mahfuz/CRISPR/2_HAP1/GIN_Analysis/FLEX/R/Pathway_Test')

GOBP.DepMap.19Q2.set1 <- CalculatePredictionAndTrueOnLibraryProfiles (data.ca, data.interaction.set1)
save(GOBP.DepMap.19Q2.set1, file = 'GOBP.DepMap.19Q2.set1.Rdata')

GOBP.DepMap.19Q2.set2 <- CalculatePredictionAndTrueOnLibraryProfiles (data.ca, data.interaction.set2)
save(GOBP.DepMap.19Q2.set2, file = 'GOBP.DepMap.19Q2.set2.Rdata')


# 5. Plot global PR Curves
# ---------------------------------------------
load('GOBP.DepMap.19Q2.set1.Rdata')
load('GOBP.DepMap.19Q2.set2.Rdata')

pred.ca <- list(out_19Q2_1 = list(true = GOBP.DepMap.19Q2.set1$true, 
                                  predicted = GOBP.DepMap.19Q2.set1$predicted))
pred.ca <- append(pred.ca, list(out_19Q2_2 = list(true = GOBP.DepMap.19Q2.set2$true, 
                                                  predicted = GOBP.DepMap.19Q2.set2$predicted)))

PlotPRSimilarity (pred.ca, subsample = TRUE, type.plot = 'log',
                  fig.title = 'DepMap 19Q2 subsets', 
                  fig.labs = c('TP', 'Precision'), legend.names = c('Set 1', 'Set 2'), 
                  legend.color = c('#de2d26', '#3182bd'), save.figure = TRUE,
                  outfile.name = 'DepMap_19Q2_GOBP', outfile.type = 'pdf')