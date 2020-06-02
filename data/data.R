#' CORUM Mammalian complexes (only human part)
#'
#' @source https://mips.helmholtz-muenchen.de/corum/#download
#' @format A data frame with columns:
#' \describe{
#'  \item{ID}{The Complex ID as indicated by 'ComplexID' column of original data}
#'  \item{Name}{Name of the complex}
#'  \item{Genes}{Genes belonging to the complex}
#'  \item{Length}{Number of genes in a complex}
#' }
"data_complex"


#' Canonical pathways (CP) data from MsigDB (includes KEGG, Reactome, BIOCARTA etc. pathways)
#'
#' @source https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#C2
#' @format A data frame with columns:
#' \describe{
#'  \item{ID}{Pathway ID: pseduo generated and doesn't hold any special meaning}
#'  \item{Name}{Name of the pathway}
#'  \item{Genes}{Genes belonging to the pathway}
#'  \item{Length}{Number of genes in a pathway}
#' }
"data_pathway"


#' GO biological processes (BP)
#'
#' @source http://geneontology.org/docs/download-ontology/
#' @format A data frame with columns:
#' \describe{
#'  \item{ID}{GO ID: an unique ID for each GO term}
#'  \item{Name}{Name of the GO term}
#'  \item{Genes}{Genes belonging to the GO term}
#'  \item{Length}{Number of genes in the GO term}
#' }
"data_GO_BP"


#' Two small subsets of the DepMap 19Q2 data (for testing)
#'
#' @source https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#C2
#' @format A data frame with genes on the row and cell lines on the column. Each cell quantified the dependency of the gene on the cell line.
"data_depmap_19Q2_test"
