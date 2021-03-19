#' Importing Histological Annotations
#'
#' ImportHistologicalAnnotations()

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
