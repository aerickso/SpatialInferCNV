#' Importing Histological Annotations
#'
#' ImportHistologicalAnnotations()

ImportOriginalSTCountData <- function(SectionName, InputCountFile) {
  input <- as.data.frame(read.delim(InputCountFile, row.names = 1))
  input <- rownames_to_column(input)
  input$rowname <- paste0(SectionName, "_", input$rowname)
  names(input)[1] <- "Barcode"
  return(input)
}
