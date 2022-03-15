#' Preparing a matrix for spatial visualization of number of genes with an inferred CNV, derived from spatial transriptomics data.
#'
#' Output_PGA_Visualization_MatrixGreyNA()
#' 
#' @param SectionName A character string for section name.
#' @param InputCNVs An input dataframe created by the function ExtractSectionWise()
#' @param BarcodesFile A single column dataframe comprised of a list of barcode coordinates in the form AxB, where A = the X coordinate, and B = the Y coordinate.
#' 
#' @return A dataframe for spatial visualization by Plot_PGA_Visualization_Matrix()
#' @examples
#' Output_PGA_Visualization_MatrixGreyNA("H2_1", H2_1_Sectionwise_CNVsGenes_Counted, L2_Barcodes)

Output_PGA_Visualization_MatrixGreyNA <- function(SectionName, InputCNVs, BarcodesFile) {
  input <- InputCNVs %>%
    mutate(section = substr(Barcode, 1, 4)) %>%
    filter(section == paste0(SectionName)) %>%
    mutate(Barcode = substr(Barcode, 6, max(nchar(Barcode))))
  PGA_Visualization_Matrix <- left_join(BarcodesFile, input)
  PGA_Visualization_Matrix$Test <- sub('.*\\_', '', PGA_Visualization_Matrix$Barcode)
  PGA_Visualization_Matrix$x <- as.numeric(sub('.*x', '', PGA_Visualization_Matrix$Test))
  PGA_Visualization_Matrix$y <- as.numeric(sub('x.*', '', PGA_Visualization_Matrix$Test))
  PGA_Visualization_Matrix <- PGA_Visualization_Matrix %>%
    select(x, y, PercentageGenomeAltered) %>%
    arrange(x, y)
  return(PGA_Visualization_Matrix)
}