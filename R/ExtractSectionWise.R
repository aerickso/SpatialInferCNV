#' ExtractSectionWise
#'
#' ExtractSectionWise()

ExtractSectionWise <- function(SectionName, CNV_Genes_Organscale_Input, AllBarcodes, Threshold) {
  output <- CNV_Genes_Organscale_Input %>%
    filter(section == paste0(SectionName)) %>%
    select(-section)
  Counted <- output %>% group_by(gene) %>% tally()
  sectionbarcodes <- AllBarcodes %>%
    filter(Histology == paste0(SectionName))
  MaxLength <- as.numeric(nrow(sectionbarcodes))
  CountPercentageThreshold <- round(Threshold * MaxLength,0)
  CountedThresholded <- Counted %>%
    filter(n > CountPercentageThreshold)
  CNV_Genes_Filtered <- inner_join(output, CountedThresholded)
  CNVs <- CNV_Genes_Filtered
  CNVsGenes_Counted <- CNVs %>% group_by(cell_group_name) %>% tally()
  names(CNVsGenes_Counted)[1] <- "Barcode"
  names(CNVsGenes_Counted)[2] <- "PercentageGenomeAltered"
  return(CNVsGenes_Counted)
}
