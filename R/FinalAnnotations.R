#' Creating A finalized annotation dataframe containing only barcodes in the count file.  
#'
#' FinalAnnotations()
#' 
#' @param InputOriginalAnnotationFile A dataframe of barcodes selected for analysis
#' @param InputCounts A joined count dataframe, of barcodes selected for analysis AND has passed QC (counts per spot >= 500 counts)
#' @return A finalized annotation dataframe containing only barcodes in the count file.
#' @examples
#' SelectingSubTreeData(my.subtrees, 4617)
#' FinalAnnotations(MergedAll, Counts_joined)

FinalAnnotations <- function(InputOriginalAnnotationFile, InputCounts) {
  input <- InputCounts
  input <- as.data.frame(input[1,])
  input <- as.data.frame(t(input))
  input <- rownames_to_column(input, var = "Barcode")
  input <- as.data.frame(input[,1])
  names(input)[1] <- "Barcode"
  input <- right_join(InputOriginalAnnotationFile, input)
  return(input)
}
