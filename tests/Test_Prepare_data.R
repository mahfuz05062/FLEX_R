## ============================ ##
## Prepare the necessary data for the package
## ============================ ##

## 1. Make two small subsets of the DepMap 19Q2 data
file.int <- '/project/chadm/Mahfuz/Database/For_Human_GI_Project/DepMap/Depmap_19Q2/Ceres_score_19Q2_with_depmap_id.txt'
data.interaction <- GetInteractionData(file.int)

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


## 2. Save data from complex standard (will be necessary to create co-annotation standards)
filename.complex <- '/project/chadm/Mahfuz/Database/Custom/Complex/v3.0/allComplexes_Human_Name_ID_Length_Genes.txt'
data_complex <- read.table(filename.complex, header=TRUE, sep="\t", quote = '', stringsAsFactors = FALSE)
save(data_complex, file = 'data_complex.rda')


## 3. Save data from pathway standard
filename.pathway <- '/project/chadm/Mahfuz/Database/Custom/Pathway/MsigDB/All/allPathways_CP_Name_ID_Genes_Length.txt'
data_pathway <- read.table(filename.pathway, header = T, sep = '\t', quote = '', stringsAsFactors = F)
save(data_pathway, file = 'data_pathway.rda')

## 4. Save data for GO BP
filename.GO.BP <- '/project/chadm/Mahfuz/Database/Custom/GO/GOIDs_TKOv3_DepMap_19Q2_Name_ID_Length_Genes.txt' 
data_GO_BP <- read.table(filename.GO.BP, header = T, sep = '\t', quote = '', stringsAsFactors = F)
save(data_GO_BP, file = 'data_GO_BP.rda')

## ============================ ##
## Miscellaneous
## ============================ ##

# Prepare DepMap data for FLEX (for the first time)
filename <- '/project/chadm/Mahfuz/Database/For_Human_GI_Project/DepMap/Depmap_18Q3/gene_effect.csv'
data.interaction <- ReformatDepMapData(filename)
save(data.interaction, file = 'Depmap_18Q3.Rdata')
