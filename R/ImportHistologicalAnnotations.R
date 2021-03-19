#' Importing Histological Annotations
#'
#' ImportHistologicalAnnotations()

ImportHistologicalAnnotations <- function(SectionName, InputAnnotationFile) {
  input <- read.csv(paste0(InputAnnotationFile))
  names(input)[2] <- "Histology"
  input <- input %>%
    mutate(Barcode = str_replace_all(Barcode, "-", "."))
  input$Barcode <- paste0(SectionName, "_", input$Barcode)
  return(input)
}
