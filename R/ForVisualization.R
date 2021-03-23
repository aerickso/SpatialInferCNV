#' Formatting Data for Visualization
#'
#' ForVisualization()

ForVisualization <- function(InputOriginalAnnotationFile, SectionName) {
  output <- InputOriginalAnnotationFile %>%
    mutate(Barcode, paste0(substr(Barcode, start = 1, stop = 4))) %>%
    rename("Section" = 3) %>%
    mutate(Barcode, paste0(substr(Barcode, start = 6, stop = 100))) %>%
    rename("Original_Barcode" = 4) %>%
    mutate(Original_Barcode = str_replace_all(Original_Barcode, "[.]", "-")) %>%
    select(Original_Barcode, Histology, Section) %>%
    filter(Section == paste0(SectionName)) %>%
    select(Original_Barcode, Histology)
  return(output)
}
