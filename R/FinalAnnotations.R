#' Importing Histological Annotations
#'
#' ImportHistologicalAnnotations()

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
