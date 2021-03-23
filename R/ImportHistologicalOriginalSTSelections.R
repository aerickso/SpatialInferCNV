#' ImportHistologicalOriginalSTSelections
#'
#' ImportHistologicalOriginalSTSelections()

ImportHistologicalOriginalSTSelections <- function(SectionName, InputAnnotationFile) {
  input <- read.delim(paste0(InputAnnotationFile), sep = "\t")
  input <- input %>% select(x, y)
  input$Barcode <- paste0(SectionName, "_",input$x, "x", input$y)
  input <- input %>% select(Barcode)
  return(input)
}
