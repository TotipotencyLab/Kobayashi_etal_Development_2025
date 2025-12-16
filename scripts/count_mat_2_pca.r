# count_mat <- dba_k27$binding[ , -(1:3)]
# dim(count_mat)
# sample_table <- dba_k27$samples %>% 
#   tibble::column_to_rownames(var="SampleID")

count_mat_2_pca <- function(count_mat, n_pc=3, as_log10=TRUE, ignore_all_zero=TRUE, label_sample=FALSE, 
                            sample_table=NULL, color_by=NULL, shape_by=NULL,
                            output_type=c("pca"), plot_theme=NULL){
  # output_type: pca, eigen, data
  debug=F
  
  output_type <- tolower(output_type)
  if(!output_type %in% c("all", "pca", "eigen", "data")){
    stop("Unknown output type. Only accept one of these")
  }
  
  if(ignore_all_zero){
    # remove genes with no expression in all samples
    count_mat <- count_mat[rowSums(count_mat == 0) != ncol(count_mat) , ]
  }
  
  if(as_log10){
    count_mat <- log10(count_mat)
    count_mat <- replace(count_mat, is.infinite(count_mat), 0)
  }
  
  # Calculate PC ----------------------------------------------------------------------------------
  mds <- cmdscale(dist(t(count_mat)), k=n_pc, eig=TRUE)  
  eig_pc <- mds$eig * 100 / sum(mds$eig)
  
  mds_df <- data.frame(mds$points)
  colnames(mds_df) <- paste0("PC", seq_along(mds_df))
  mds_df$sample <- rownames(mds_df)
  
  # Combine PC data with sample_table, and figure whether we want to do color_by and shape_by when plot
  pc_df <- mds_df
  do_color_by=FALSE; do_shape_by=FALSE # Set default
  if(inherits(sample_table, "data.frame")){
    if(!is.null(rownames(sample_table))){
      if(debug){
        cat("pc_df:\n----------------\n")
        print(pc_df)
        
        cat("\nsample_table:\n----------------\n")
        print(sample_table)
      }
      
      if(identical(sort(rownames(sample_table)), sort(rownames(pc_df)))){
        req_col <- c(color_by, shape_by)
        if(!is.null(req_col)){
          if(all(!req_col %in% colnames(sample_table))){
            warning("None of color_by or shape_by column exist in sample_table")
          }
          found_req_col <- req_col[req_col %in% colnames(sample_table)]
          pc_df <- cbind(pc_df, sample_table[rownames(pc_df) , found_req_col, drop=FALSE])
          # Determine how the point should be made
          if(!is.null(color_by)) if(color_by %in% colnames(pc_df)) do_color_by=TRUE
          if(!is.null(shape_by)) if(shape_by %in% colnames(pc_df)) do_shape_by=TRUE
        }else{
          warning("The parameter color_by or shape_by must be provided.")
        }
      }else{
        warning("Unmatched rownames - sample_table must have rowname that matches colnames of count_mat")
      }
    }else{
      warning("sample_table doesn't have rownames - sample_table must have rowname that matches colnames of count_mat")
    }
  }
  if(debug) print(head(pc_df))
  
  if(output_type=="data"){
    out_list <- list(pca = pc_df, eigen = eig_pc)
    return(out_list)
  }
  
  # Making PCA plot
  if(is.null(plot_theme)){
    plot_theme <- ggplot2::theme_bw() +
      ggplot2::theme(text = ggplot2::element_text(size = 15),
                     axis.line = ggplot2::element_line(colour = "black"), 
                     legend.box="vertical", legend.direction="horizontal", legend.position="bottom")
  }
  
  # PCA -------------------------------------------------------------------------------------------
  if(output_type %in% c("pca", "all")){
    stt = ifelse(as_log10, yes="log10-transformed value", no="")
    x_lab <- paste0("PC1 (", round(eig_pc[1], 2), "%)")
    y_lab <- paste0("PC2 (", round(eig_pc[2], 2), "%)")
    
    p_pca <- ggplot(pc_df, aes(x=PC1, y=PC2)) + 
      labs(title = "PCA", subtitle=stt, x=x_lab, y=y_lab) +
      plot_theme
    
    point_size=5; point_alpha=0.8
    if(do_color_by & do_shape_by){
      if(debug) cat("col & shape\n")
      p_pca <- p_pca + geom_point(aes(color=.data[[color_by]], shape=.data[[shape_by]]), size=point_size, alpha=point_alpha)
    }else if(do_color_by){
      if(debug) cat("col\n")
      p_pca <- p_pca + geom_point(aes(color=.data[[color_by]]), size=point_size, alpha=point_alpha)
    }else if(do_shape_by){
      if(debug) cat("shape\n")
      p_pca <- p_pca + geom_point(aes(shape=.data[[shape_by]]), size=point_size, alpha=point_alpha)
    }else{
      if(debug) cat("plain\n")
      p_pca <- p_pca + geom_point(size=point_size, alpha=point_alpha)
    }
    
    # Add sample label
    if(label_sample){
      # Prepare data
      if("sample" %in% colnames(pc_df)){
        warning("Column 'sample' will be replaced by the rownames of the pca_data")
        pc_df <- pc_df[ , -which(colnames(pc_df) == "sample"), drop=FALSE]
        pc_df$sample <- rownames(pc_df)
      }
      # Configure which function to use for the plot
      label_func <- geom_point
      if("ggrepel" %in% rownames(installed.packages())){
        label_func <- ggrepel::geom_text_repel
      }
      # Add label
      if(do_color_by){
        p_pca <- p_pca + label_func(aes(label=sample, color=.data[[color_by]]), max.overlaps=Inf)
      }else{
        p_pca <- p_pca + label_func(aes(label=sample), max.overlaps=Inf)
      }
    }
  } # end making PCA
  
  
  # Eigen -----------------------------------------------------------------------------------------
  if(output_type %in% c("eigen", "all")){
    eigen_df <- data.frame(pc = seq_along(eig_pc), percent = eig_pc)
    p_eigen <- ggplot(eigen_df, aes(x=pc, y=percent)) +
      geom_col() +
      labs(title = "Percent variance explained by each PC", x="PC", y="% Explained") +
      scale_x_continuous(breaks=seq_along(eig_pc)) +
      plot_theme
  }
  
  # Return output ---------------------------------------------------------------------------------
  if(output_type == "all"){
    out_list <- list(
      pca = p_pca,
      eigen = p_eigen,
      pca_data = pc_df,
      eigen_vec = eigen_df
    )
    return(out_list)
  }else if(output_type == "pca"){
    return(p_pca)
  }else if(ouput_type == "eigen"){
    return(p_eigen)
  }else{
    stop("Unrecognized output type")
  }
}

