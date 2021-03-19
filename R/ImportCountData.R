#' Importing Histological Annotations
#'
#' ImportHistologicalAnnotations()

ImportCountData <- function(SectionName, InputCountFile) {
  input <- Read10X_h5(InputCountFile, use.names = FALSE)
  input <- as.matrix(input)
  input <- as.data.frame(t(input))
  input <- rownames_to_column(input, "Barcode")
  input$Barcode <- paste0(SectionName, "_", input$Barcode)
  input$Barcode <- gsub("\\-", "\\.", input$Barcode)
  return(input)
}
