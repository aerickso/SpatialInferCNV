#' Importing histological annotations of Visium barcodes and appending a section name to the barcodes.
#'
#' ImportHistologicalAnnotations()
#' 
#' @param SectionName A character string for section name.
#' @param InputAnnotationFile A file path to a .tsv file
#' @return A dataframe of barcodes with appended section names
#' @examples
#' ImportHistologicalAnnotations("H1_2", "./H1_2_Final_Consensus_Annotations.csv")

ImportHistologicalAnnotations <- function(SectionName, InputAnnotationFile) {
  input <- read.csv(paste0(InputAnnotationFile))
  names(input)[2] <- "Histology"
  input <- input %>%
    mutate(Barcode = str_replace_all(Barcode, "-", "."))
  input$Barcode <- paste0(SectionName, "_", input$Barcode)
  return(input)
}
