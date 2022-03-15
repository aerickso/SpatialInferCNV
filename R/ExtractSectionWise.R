#' Obtaining a thresholded dataframe as part of spatial visualization of spatial transcriptomics data.
#'
#' ExtractSectionWise()
#' 
#' @param SectionName A character string for section name.
#' @param CNV_Genes_Organscale_Input A dataframe, mirroring the structure of infercnv::run output file 17_HMM_predHMMi6.hmm_mode-cells.pred_cnv_genes.dat
#' @param AllBarcodes A dataframe of barcodes and annotations.
#' @param Threshold A numerical value for sectionwise thresholding of the number of genes to pass: integer values from 0-100.
#' 
#' @return A dataframe of ST counts, that have passed QC and are selected.
#' 
#' @examples
#' ExtractSectionWise("H2_1", CNV_Genes_Filtered, AllBarcodes, 0.45)
 
 
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
