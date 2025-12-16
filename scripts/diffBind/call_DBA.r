call_DBA <- function(dds_res_df, log2FoldChange_cutoff=log2(2), padj_cutoff=0.05, 
                     overwrite_result_column = FALSE, result_colname="DBA", 
                     gain = "Gain", lost = "Lost", no_change="No_Change") {
  require(dplyr)
  # Assuming the feature identifier is in the rownames
  # This should result to data.frame with 1 column (`result_colname`) that annotate differential expressed genes
  
  # Checking inputs -------------------------------------------------------------------------------
  # Check Input length
  input_length <- sapply(list(result_colname=result_colname, gain=gain, lost=lost, no_change=no_change), length)
  if(any(input_length != 1)){
    wrong_inputs <- input_length[input_length != 1]
    msg <- "Input(s) has incorrect length. Only accept the input of length 1 for these options: \n"
    for(i in seq_along(wrong_inputs)){
      msg <- c(msg, "  ", names(wrong_inputs)[i], " - received input length: ", wrong_inputs[i], "\n")
    }
    stop(msg)
  }
  
  # Check requested DEG labels
  if(any(duplicated(c(gain, lost, no_change)))){
    label_vec <- c(gain=gain, lost=lost, no_change=no_change)
    dup_list <- split(names(label_vec), label_vec)
    dup_list <- dup_list[sapply(dup_list, length) > 1]
    msg <- "Duplicated DEG labels is not allowed\n"
    for(i in seq_along(dup_list)){
      msg <- c(msg, "[ ", paste0(dup_list[[i]], collapse=", "), " ] is labeled as '", names(dup_list)[i], "'\n")
    }
    stop(msg)
  }
  
  # Check requested result_colname
  if(result_colname %in% colnames(dds_res_df)){
    msg <- c("The column ", result_colname, " is already exist in the column of input data frame.")
    if(overwrite_result_column){
      msg <- c(msg, "\nOverwrite the results over this column")
      warning(msg)
      
      # make sure the existing column is gone
      dds_res_df <- dds_res_df[ , colnames(dds_res_df) != result_colname, drop=FALSE]
    }else{
      stop(msg)
    }
  }
  
  # Actual DEG annotation -------------------------------------------------------------------------
  deg_df <- dds_res_df %>% 
    dplyr::transmute(
      sig = !is.na(padj) & (padj < padj_cutoff) & (abs(log2FoldChange) > abs(log2FoldChange_cutoff)),
      DBA = dplyr::case_when(!sig ~ no_change,
                             (log2FoldChange > 0) ~ gain,
                             (log2FoldChange < 0) ~ lost, 
                             .default=no_change)
    ) %>% 
    dplyr::select(-sig) %>% 
    dplyr::rename(!!sym(result_colname) := DBA)
  
  dds_res_df <- cbind(dds_res_df, deg_df)
  return(dds_res_df)
}
