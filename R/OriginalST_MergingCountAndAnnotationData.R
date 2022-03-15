#' Merging spatial transcriptomics, 1k array count files and barcodes, an apply a QC metric to only select
#' ST spots with >=500 total unique molecular identifiers.
#'
#' OriginalST_MergingCountAndAnnotationData()
#' 
#' @param InputAnnotationFile An annotation file created by ImportHistologicalOriginalSTSelections()
#' @param InputCountFile A ST count file created by ImportOriginalSTCountData()
#' @return A dataframe of ST counts, that have passed QC and are selected.
#' @examples
#' OriginalST_MergingCountAndAnnotationData(Barcodes_H2_1, Counts_H2.1)

OriginalST_MergingCountAndAnnotationData <- function(InputAnnotationFile, InputCountFile) {
  formerge <- select(InputAnnotationFile, -Histology)
  MergedAnnotationsandCounts <- inner_join(formerge, InputCountFile)
  MergedAnnotationsandCounts <- remove_rownames(MergedAnnotationsandCounts)
  MergedAnnotationsandCounts <- column_to_rownames(MergedAnnotationsandCounts, "Barcode")
  MergedAnnotationsandCounts$Total <- rowSums(MergedAnnotationsandCounts)
  MergedAnnotationsandCounts <- MergedAnnotationsandCounts %>% filter(Total >= 500)
  MergedAnnotationsandCounts <- select(MergedAnnotationsandCounts, -Total)
  MergedAnnotationsandCounts <- as.data.frame(t(MergedAnnotationsandCounts))
  if(length(MergedAnnotationsandCounts) == 1){
    MergedAnnotationsandCounts <- tibble::rownames_to_column(MergedAnnotationsandCounts, "Genes")
    return(MergedAnnotationsandCounts)
  } else {
    MergedAnnotationsandCounts <- MergedAnnotationsandCounts[,colSums(is.na(MergedAnnotationsandCounts))<nrow(MergedAnnotationsandCounts)]
    MergedAnnotationsandCounts <- tibble::rownames_to_column(MergedAnnotationsandCounts, "Genes")
    return(MergedAnnotationsandCounts)
  }
}
