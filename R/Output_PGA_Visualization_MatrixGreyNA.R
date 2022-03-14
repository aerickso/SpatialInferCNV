#' Output_PGA_Visualization_MatrixGreyNA
#'
#' Output_PGA_Visualization_MatrixGreyNA()

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