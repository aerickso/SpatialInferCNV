#' Importing Histological Annotations
#'
#' ImportHistologicalAnnotations()

MappingGeneNames <- function(MappingFile, CountsJoinedObject) {
  GenesForMapping <- MappingFile %>% select(ENSMBLID, chr, start, stop)
  GenesInSample <- CountsJoinedObject %>% select(Genes)
  GenesInSamplevsOrdering <- merge(GenesInSample, GenesForMapping)
  dedup_GenesInSamplevsOrdering <- GenesInSamplevsOrdering[!duplicated(GenesInSamplevsOrdering$Genes), ]
  dedup_GenesInSamplevsOrdering$chromorder <- gsub("chr","",dedup_GenesInSamplevsOrdering$chr)
  table(dedup_GenesInSamplevsOrdering$chromorder)
  dedup_GenesInSamplevsOrdering$chromorder <- as.numeric(ifelse(dedup_GenesInSamplevsOrdering$chromorder == "X", 23,
                                                                ifelse(dedup_GenesInSamplevsOrdering$chromorder == "Y", 24,      dedup_GenesInSamplevsOrdering$chromorder)))
  dedup_GenesInSamplevsOrdering <- dedup_GenesInSamplevsOrdering[order(dedup_GenesInSamplevsOrdering$chromorder),]
  dedup_GenesInSamplevsOrdering <- dedup_GenesInSamplevsOrdering[,1:4]
  return(dedup_GenesInSamplevsOrdering)
}
