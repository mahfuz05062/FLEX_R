## ============================ ##
## Prepare the necessary data for the package
## ============================ ##

## 1. Make two small subsets of the DepMap 19Q2 data
# file.int <- '/home/mahfuz/Desktop/CRISPR/Datasets/DepMap_19Q2/Ceres_score_19Q2_with_depmap_id.txt'
file.int <- '/project/chadm/Mahfuz/Database/For_Human_GI_Project/DepMap/Depmap_19Q2/Ceres_score_19Q2_with_depmap_id.txt'
data.interaction <- GetInteractionData(file.int)

# tmp <- apply(data.interaction, 1, sd)

# Set 1
data.interaction.set1 <- data.interaction[, 1:50]
if (is.unsorted(rownames(data.interaction.set1))){ # Sort by gene names
  ind <- order(rownames(data.interaction.set1))
  data.interaction.set1 <- data.interaction.set1[ind,]
}

# Set 2
data.interaction.set2 <- data.interaction[, 51:100]
if (is.unsorted(rownames(data.interaction.set2))){
  ind <- order(rownames(data.interaction.set2))
  data.interaction.set2 <- data.interaction.set2[ind,]
}

save(data.interaction.set1, data.interaction.set2, file = 'data_depmap_19Q2_test.rda')

# scp rahma118@mcclintock.cs.umn.edu:/project/chadm/Mahfuz/CRISPR/2_HAP1/GIN_Analysis/FLEX/R/FLEX/data/*.rda /home/mahfuz/Desktop/CRISPR/FLEX/R/FLEX/data/

## 2. Save data from complex standard (will be necessary to create co-annotation standards)
filename.complex <- '/project/chadm/Mahfuz/Database/Custom/Complex/v3.0/allComplexes_Human_Name_ID_Length_Genes.txt'
data_complex <- read.table(filename.complex, header=TRUE, sep="\t", quote = '', stringsAsFactors = FALSE)
save(data_complex, file = 'data_complex.rda')


## 3. Save data from pathway standard
setwd('/home/mahfuz/Desktop/CRISPR/FLEX/Matlab')
filename.pathway <- '/project/chadm/Mahfuz/Database/Custom/Pathway/MsigDB/All/allPathways_CP_Name_ID_Genes_Length.txt'
data_pathway <- read.table(filename.pathway, header = T, sep = '\t', quote = '', stringsAsFactors = F)
save(data_pathway, file = 'data_pathway.rda')

# scp /home/mahfuz/Desktop/CRISPR/FLEX/R/FLEX/tests/*.R rahma118@mcclintock.cs.umn.edu:/project/chadm/Mahfuz/CRISPR/2_HAP1/GIN_Analysis/FLEX/R/FLEX/tests/

## 4. Save data for GO BP
filename.GO.BP <- '/project/chadm/Mahfuz/Database/Custom/GO/GOIDs_TKOv3_DepMap_19Q2_Name_ID_Length_Genes.txt' 
data_GO_BP <- read.table(filename.GO.BP, header = T, sep = '\t', quote = '', stringsAsFactors = F)
save(data_GO_BP, file = 'data_GO_BP.rda')

# scp rahma118@mcclintock.cs.umn.edu:/project/chadm/Mahfuz/CRISPR/2_HAP1/GIN_Analysis/FLEX/R/FLEX/data/data_GO_BP.rda /home/mahfuz/Desktop/CRISPR/FLEX/R/FLEX/data/

## ============================ ##
## Miscellaneous
## ============================ ##

# Prepare DepMap data for FLEX (for the first time)
filename <- '/project/chadm/Mahfuz/Database/For_Human_GI_Project/DepMap/Depmap_18Q3/gene_effect.csv'
data.interaction <- ReformatDepMapData(filename)
save(data.interaction, file = 'Depmap_18Q3.Rdata')