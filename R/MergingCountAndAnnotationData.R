#' Merging Visium spatial transciptomics count and annotation data, as well as applying a QC filter to only include spots with >= 500 counts
#'
#' MergingCountAndAnnotationData()
#' 
#' @param SectionName A character string for section name.
#' @param InputAnnotationFile An annotation file containing all barcodes to be used in the analysis (bound dataframe of one or more outputs from ImportHistologicalAnnotations())
#' @param InputCountFile A dataframe of Visium count data (output from ImportCountData())
#' @return A dataframe of barcodes with appended section names that have passed QC
#' @examples
#' MergingCountAndAnnotationData("H2_1",MergedAll, H2_1_ENSBMLID_Counts)

MergingCountAndAnnotationData <- function(SectionName, InputAnnotationFile, InputCountFile) {
  formerge <- select(InputAnnotationFile, -Histology)
  MergedAnnotationsandCounts <- inner_join(formerge, InputCountFile)
  MergedAnnotationsandCounts <- remove_rownames(MergedAnnotationsandCounts)
  MergedAnnotationsandCounts <- column_to_rownames(MergedAnnotationsandCounts, "Barcode")
  MergedAnnotationsandCounts$Total <- rowSums(MergedAnnotationsandCounts)
  MergedAnnotationsandCounts <- MergedAnnotationsandCounts %>% filter(Total >= 500)
  MergedAnnotationsandCounts <- select(MergedAnnotationsandCounts, -Total)
  MergedAnnotationsandCounts <- as.data.frame(t(MergedAnnotationsandCounts))
  MergedAnnotationsandCounts <- MergedAnnotationsandCounts[,colSums(is.na(MergedAnnotationsandCounts))<nrow(MergedAnnotationsandCounts)]
  MergedAnnotationsandCounts <- tibble::rownames_to_column(MergedAnnotationsandCounts, "Genes")
  return(MergedAnnotationsandCounts)
}
