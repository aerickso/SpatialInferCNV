#' Importing Histological Annotations
#'
#' ImportHistologicalAnnotations()

ImportGEFAnnotations <- function(SectionName, InputAnnotationFile) {
  input <- read.csv(paste0(InputAnnotationFile))
  input$Barcode <- paste0(SectionName, "_", input$Barcode)
  input$Barcode <- gsub("\\-", "\\.", input$Barcode)
  return(input)
}
