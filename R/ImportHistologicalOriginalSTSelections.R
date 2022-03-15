#' Importing spatial transcriptomics, 1k array selected spot file data and append section names to the barcodes.
#'
#' ImportHistologicalOriginalSTSelections()
#' 
#' @param SectionName A character string for section name.
#' @param InputAnnotationFile A file path to a .tsv file
#' @return A dataframe of barcodes with appended section names
#' @examples
#' ImportHistologicalOriginalSTSelections("H2_1", "./Patient 1/1k_arrays/H2_1/spot_data-selection-180903_L11_CN63_D1_P_H2.1_CY3_EB_aligned.tsv")

ImportHistologicalOriginalSTSelections <- function(SectionName, InputAnnotationFile) {
  input <- read.delim(paste0(InputAnnotationFile), sep = "\t")
  input <- input %>% select(x, y)
  input$Barcode <- paste0(SectionName, "_",input$x, "x", input$y)
  input <- input %>% select(Barcode)
  return(input)
}
