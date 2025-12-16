dba_to_dds <- function(dba){
  # Get row and column data
  row_gr <- dba_get_merged_peak(dba, add_peak_call=TRUE)
  names(row_gr) <- paste0("peak_", seq_along(row_gr))
  row_df <- as.data.frame(row_gr)
  col_df <- tibble::column_to_rownames(dba$samples, var="SampleID")
  
  # export to dds object for differential binding analysis (for < 3 reps)
  m <- matrix(NA, nrow=length(row_gr), ncol=nrow(col_df),
              dimnames=list(names(row_gr), rownames(col_df)))
  for(i in seq_along(dba$peaks)){
    # dba$peaks is unnamed, assume it has the same name as sample
    m[ , colnames(dba$class)[i]] <- dba$peaks[[i]]$Reads
  }
  
  # Construct dds object
  dds <- DESeqDataSetFromMatrix(countData=m, colData=col_df, rowData=row_gr, design=as.formula("~Condition"))
  
  # inherit size factor
  if("norm" %in% names(dba)){
    metadata(dds)$norm <- dba$norm
    sizeFactors(dds) <- dba$norm$DESeq2$norm.facs
  }
  
  return(dds)
}
