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


## 5. Download data for Functional network (GIANT)
#  See 'MakeFuncNetFromGIANT' in GoldStandards.R
# https://www.programmableweb.com/news/how-to-access-any-restful-api-using-r-language/how-to/2017/07/21

# install.packages("httr")
# install.packages("jsonlite")
# require("httr")
# require("jsonlite")

# Getting the data
# base <- 'https://hb.flatironinstitute.org/api/datasets' # Not this one (https://hb.flatironinstitute.org/data)

# base <- 'https:/hb.flatironinstitute.org/api/integrations/'# We want networks
# get_prices <- GET(base)
# get_prices_text <- content(get_prices, "text")
# get_prices_json <- fromJSON(get_prices_text, flatten = TRUE)
# get_prices_df <- as.data.frame(get_prices_json)
# out <- get_prices_df[which(get_prices_df$slug == 'global'), ]

# How to get only the top edge file?
# file_location <- '/project/chadm/Mahfuz/Database/Greene_et_al/Func_net/global_top.gz'
# download.file(url='https://s3-us-west-2.amazonaws.com/humanbase/networks/global_top.gz', destfile = file_location, method='curl') # This is the direct link (but I would rather love to use the api to get this link!)


## ============================ ##
## Miscellaneous
## ============================ ##

# Prepare DepMap data for FLEX (for the first time)
filename <- '/project/chadm/Mahfuz/Database/For_Human_GI_Project/DepMap/Depmap_18Q3/gene_effect.csv'
data.interaction <- ReformatDepMapData(filename)
save(data.interaction, file = 'Depmap_18Q3.Rdata')
