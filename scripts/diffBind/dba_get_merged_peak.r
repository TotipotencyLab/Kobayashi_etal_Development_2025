dba_get_merged_peak <- function(dba, add_peak_call=FALSE){
  merged_peaks <- as.data.frame(dba$merged)
  colnames(merged_peaks) <- tolower(colnames(merged_peaks))
  merged_peaks$chr <- dba$chrmap[dba$merged[ , 1]]
  
  if(add_peak_call){
    merged_peaks <- cbind(merged_peaks, dba$called>0)
  }
  
  merged_peaks_gr <- makeGRangesFromDataFrame(merged_peaks, keep.extra.columns=TRUE)
  return(merged_peaks_gr)
}
