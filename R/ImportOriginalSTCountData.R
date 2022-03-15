#' Importing spatial transcriptomics, 1k array count data and append section names to the barcodes.
#'
#' ImportOriginalSTCountData()
#'
#' @param SectionName A character string for section name.
#' @param InputCountFile A file path to a .tsv file
#' @return A dataframe of count data, having barcodes with appended section names
#' @examples
#' ImportOriginalSTCountData("H2_1", "./Patient 1/1k_arrays/H2_1/180903_L11_CN63_D1_H2.1_EB_stdata.tsv")

ImportOriginalSTCountData <- function(SectionName, InputCountFile) {
  input <- as.data.frame(read.delim(InputCountFile, row.names = 1))
  input <- rownames_to_column(input)
  input$rowname <- paste0(SectionName, "_", input$rowname)
  names(input)[1] <- "Barcode"
  return(input)
}
