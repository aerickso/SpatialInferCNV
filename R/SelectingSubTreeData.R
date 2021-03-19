#' Selecting Subtree Data for Node Selection
#'
#' SelectingSubTreeData()

SelectingSubTreeData <- function(SubtreeObject, NodeOfInterest) {
  tree_node <- SubtreeObject[[NodeOfInterest]]
  output <- tree_node$tip.label
  output <- as.data.frame(output)
  output <- output %>%
    mutate(Node = paste0("Node_", NodeOfInterest))
  names(output)[1] <- "Barcode"
  return(output)
}
