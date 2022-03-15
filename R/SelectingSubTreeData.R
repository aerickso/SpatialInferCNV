#' Selecting Subtree Data for Node Selection: this selects a number of barcoded spots from a inferCNV dendrogram object for further analysis.
#'
#' SelectingSubTreeData()
#' 
#' @param SubtreeObject A dendrogram, phylo object created by subtrees(as.phylo([dendogram.txt]))
#' @param NodeOfInterest A numerical integer corresponding to a phylogram/dendogram node of interest
#' @return A specific subtree node
#' @examples
#' SelectingSubTreeData(my.subtrees, 4617)

SelectingSubTreeData <- function(SubtreeObject, NodeOfInterest) {
  tree_node <- SubtreeObject[[NodeOfInterest]]
  output <- tree_node$tip.label
  output <- as.data.frame(output)
  output <- output %>%
    mutate(Node = paste0("Node_", NodeOfInterest))
  names(output)[1] <- "Barcode"
  return(output)
}
