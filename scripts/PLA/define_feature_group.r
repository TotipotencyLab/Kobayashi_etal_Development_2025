# Identify a group of ROIs across z-stacks that potentially represent part of the cellular feature/compartment.
#' Define Feature Groups from ROI Data
#'
#' @description
#' This function takes regions of interest (ROI) defined from FIJI, then define cellular feature based on how they overlap across the z-stacks
#' The function was optimized to detect nucleus of embryos and may not work well with some other type of feature
#'
#' @param roi_df A data frame containing ROI information with spatial coordinates and z-stack positions.
#'
#' [ROI filtering parameters]
#' @param pre_roi_filter_colname Character. A column name of `roi_df` containing logical value for filtering ROI.
#' @param roi_area_range Numerical vector (length 2). A range of minimum and maximum area of ROI that should be included in the analysis. 
#'   Default: c(0, Inf), that is, no RI filtering by area.
#' @param roi_regex Character or NULL. Regular expression of ROI name to be included for the analysis. Default: NULL, no ROI filtering by name
#'
#' [ROI overlap filtering parameters]
#' @param min_intersect_ratio Numeric, between 0 and 1. Minimum intersection ratio between polygons required for consider them as part of the same feature.
#'   The intersection ratio is defined as a proportion of overlapping area and the area of smaller ROIs. 
#'   Default: 0, no filter of ROI interaction by overlapping area.
#' 
#' [Feature filtering parameters]
#' @param max_z_dist Numeric. Maximum z-distance allowed between the next overlaping ROIs to be considered part of the same feature group.
#' @param min_z_span Numeric. Minimum z-span (range across z-stacks) required for a feature group to be retained.
#' @param min_avg_area Numeric. Minimum average area threshold for filtering feature groups.
#' 
#' [Output parameters]
#' @param feature_prefix Character. Prefix of the newly defined feature
#' @param invalid_feature_prefix Character. Prefix of the newly defined feature that are not considered as valid features.
#' @param fail_ROI_feature_prefix Character. Prefix of the ROI that didn't used for defining features
#' 
#' @details
#' The function makes two key assumptions:
#' \itemize{
#'   \item Each ROI is located on only one z-stack
#'   \item No overlap exists between polygons within the same z-stack (TODO: consider addressing this in the future)
#' }
#'
#' The function depends on the following helper functions from
#' \code{FnGroup_roi_2_polygons.r}:
#' \itemize{
#'   \item \code{roi_extract_xy_coord}: Extracts XY coordinates from ROI data
#'   \item \code{roi_2_polygons}: Converts ROI data to polygon objects
#'   \item \code{polygonize_roi_df}: Polygonizes ROI dataframe
#' }
#'
#' @return A processed data structure containing defined feature groups after
#'   applying all filtering criteria.
#'
#' @seealso \code{\link{roi_extract_xy_coord}}, \code{\link{roi_2_polygons}},
#'   \code{\link{polygonize_roi_df}}



define_feature_group <- function(roi_df, 
                                 # ROI filtering
                                 pre_roi_filter_colname = "include", roi_area_range = c(0, Inf), roi_regex = NULL,
                                 # ROI overlapping filter
                                 min_intersect_ratio=0.0,
                                 # Feature filtering
                                 max_z_dist=1, min_z_span=5, min_avg_area=NULL, 
                                 # Output control
                                 feature_prefix = "feature_", invalid_feature_prefix = "invalid_feature_", fail_ROI_feature_prefix = "failed_ROI_",
                                 verbose = FALSE){
  
  # Input processing ------------------------------------------------------------------------------
  reject_input <- FALSE
  if(all(c("x", "y", "z") %in% colnames(roi_df))){
    # Collapse x,y,z coordinate in each ROI into polygon objects (sf package)
    roi_pg_df <- polygonize_roi_df(roi_df)
    
  }else if(("geometry" %in% colnames(roi_df))){
    if(("sfc_POLYGON" %in% class(roi_df$geometry))){
      # i.e., the input already has polygon
      roi_pg_df <- roi_df
    }else{
      reject_input <- TRUE
    }
  }else{
    reject_input <- TRUE
  }
  
  if(reject_input){stop("Incorrect type of input")}
  
  
  roi_pg_df$area <- st_area(roi_pg_df$geometry) # Assign area
  
  # Making a reference vector with ROI ID as a key for later on when preparing the output
  roi_z_map <- dplyr::pull(roi_pg_df, z, name=roi)
  roi_area_map <- dplyr::pull(roi_pg_df, area, name=roi)
  
  # Filtering out ROI -----------------------------------------------------------------------------
  
  ## By pre-determined excluded ROI ----
  roi_fail_df <- tibble()
  if(pre_roi_filter_colname %in% colnames(roi_df)){
    # i.e. user pre-filter the ROI before using with this function
    use_roi_flag <- roi_df[[pre_roi_filter_colname]]
    if(sum(!use_roi_flag)>0){
      valid_roi <- roi_df$roi[use_roi_flag]
      # The unused ROI will be excluded from the analysis but will merge back in the output
      roi_fail_df <-  dplyr::filter(roi_pg_df, !(roi %in% valid_roi)) %>% 
        mutate(feature_id = paste0(fail_ROI_feature_prefix, "excluded")) %>% 
        rbind(roi_fail_df, .)
      
      roi_pg_df <- dplyr::filter(roi_pg_df, (roi %in% valid_roi))
      
      if(nrow(roi_pg_df) == 0){
        warning("All ROIs were filtered out (include flag)")
        return(roi_fail_df)
      }
    }
  }
  
  ## By out by ROI name ----
  if(!is.null(roi_regex)){
    # Expect roi_regex to have a length of 1
    valid_roi <- roi_pg_df %>% 
      dplyr::filter(!is.na(roi)) %>% 
      dplyr::filter(str_detect(roi, roi_regex)) %>% 
      dplyr::pull(roi) %>% unique()
    
    roi_fail_df <- dplyr::filter(roi_pg_df, !(roi %in% valid_roi)) %>% 
      mutate(feature_id = paste0(fail_ROI_feature_prefix, "name")) %>% 
      rbind(roi_fail_df, .)
    
    roi_pg_df <- dplyr::filter(roi_pg_df, (roi %in% valid_roi))
    
    if(nrow(roi_pg_df) == 0){
      warning("All ROIs were filtered out (ROI name filtering)")
      return(roi_fail_df)
    }
  }
  
  # By by area 
  if(!is.null(roi_area_range)){
    # Expect roi_area_range to have a length of 2
    valid_roi <- roi_pg_df %>% 
      dplyr::filter(area >= min(roi_area_range, na.rm=TRUE),
                    area <= max(roi_area_range, na.rm=TRUE),
                    !is.na(roi)) %>% 
      dplyr::pull(roi) %>% unique()
    
    roi_fail_df <- dplyr::filter(roi_pg_df, !(roi %in% valid_roi)) %>% 
      mutate(feature_id = paste0(fail_ROI_feature_prefix, "area")) %>% 
      rbind(roi_fail_df, .)
    
    roi_pg_df <- dplyr::filter(roi_pg_df, (roi %in% valid_roi))
    
    if(nrow(roi_pg_df) == 0){
      warning("All ROIs were filtered out (ROI Area)")
      return(roi_fail_df)
    }
  }
  
  
  # for internal tracking of ROI at each stage of the script
  roi_set <- list(
    start = unique(dplyr::pull(dplyr::filter(roi_df, !is.na(roi)), roi)),
    input = unique(dplyr::pull(dplyr::filter(roi_pg_df, !is.na(roi)), roi))
  )
  
  # Find overlap ROI ------------------------------------------------------------------------------
  ovl_pair_df <- find_ROI_z_intersect(roi_pg_df, max_z_dist=max_z_dist, min_intersect_ratio=min_intersect_ratio, verbose=verbose)
  ovl_pair_df <- dplyr::filter(ovl_pair_df, !is.na(roi_1), !is.na(roi_2)) # Just in case
  roi_set$overlap <- unique(c(ovl_pair_df$roi_1, ovl_pair_df$roi_2))
  
  # Which ROI got removed out from this step
  roi_overlap_fail <- roi_set$input[!(roi_set$input %in% roi_set$overlap)]
  if(length(roi_overlap_fail) > 0){
    roi_fail_df <- dplyr::filter(roi_pg_df, (roi %in% roi_overlap_fail)) %>% 
      mutate(feature_id = paste0(fail_ROI_feature_prefix, "overlap")) %>% 
      rbind(roi_fail_df, .)
    
    roi_pg_df <- dplyr::filter(roi_pg_df, !(roi %in% roi_overlap_fail))
  }
  
  # Grouping features (consider making a function) ------------------------------------------------
  
  # Borrow the terminology from the graph network field
  #   node: each ROI
  #   edge: overlap between ROI
  
  # This is similar to roi_pg_df, but without the geometry field
  # Doing it this way will include ROIs that are filtering out when finding z-overlap
  # roi_pg_df %>% 
  #   dplyr::select(roi_id = roi, z, area) %>% 
  #   mutate(feature_group = NA)
  
  roi_node_df <- data.frame(roi_id = unique(c(ovl_pair_df$roi_1, ovl_pair_df$roi_2))) %>%
    as_tibble() %>% 
    mutate(z = roi_z_map[roi_id], 
           feature_group=NA, # place holder
           area = roi_area_map[roi_id])
  
  roi_edge_df <- ovl_pair_df %>%
    dplyr::select(roi_1, roi_2)
  roi_edge_mat <- as.matrix(roi_edge_df)
  
  # Assigning feature groups
  feature_count <- 0
  for(i in 1:nrow(roi_edge_mat)){
    cur_idx <- which(roi_node_df$roi_id %in% roi_edge_mat[i, ])
    # Check if any of these two already has group number assign to it
    cur_nuc_num <- unique(roi_node_df$feature_group[cur_idx]) %>% 
      subset(., !is.na(.))
    
    # Assign nucleus ID
    if(length(cur_nuc_num)==0){
      # i.e., new nucleus found!
      feature_count <- feature_count +1
      roi_node_df$feature_group[cur_idx] <- feature_count
      
    }else if(length(cur_nuc_num)==1){
      # Adding new ROI to the existing group
      roi_node_df$feature_group[cur_idx] <- cur_nuc_num
      
    }else if(length(cur_nuc_num)==2){
      # Joining the two assigned nucleus together
      cur_idx <- which(roi_node_df$feature_group %in% cur_nuc_num)
      roi_node_df$feature_group[cur_idx] <- min(cur_nuc_num)
    }
  }
  
  # Filtering feature -----------------------------------------------------------------------------
  # By number of z-span
  # Calculate number of z-span per nucleus
  z_span_df <- roi_node_df %>% 
    dplyr::filter(!is.na(feature_group)) %>% 
    dplyr::select(feature_group, z) %>% 
    unique() %>% 
    {table(.$feature_group)} %>% 
    as.data.frame() %>% as_tibble() %>% 
    `colnames<-`(c("feature_group", "n_z_span")) %>% 
    mutate(feature_group = as.double(feature_group))
  
  valid_feature_group <- z_span_df %>% 
    dplyr::filter(n_z_span >= min_z_span) %>% 
    pull(feature_group)
  
  # Filter by average area
  if(!is.null(min_avg_area)){
    nuc_stats_df <- roi_node_df %>% 
      group_by(feature_group) %>% 
      reframe(mean_area = mean(area))
    
    valid_feature_group_byArea <- nuc_stats_df %>% 
      dplyr::filter(mean_area >= min_avg_area) %>% 
      pull(feature_group) 
    
    # Filtering the list
    valid_feature_group <-  base::intersect(valid_feature_group, valid_feature_group_byArea)
  }
  
  ## Assign feature_id --------------------------------------------------------------------------------
  # valid nuc
  valid_feature_id_map <- paste0(feature_prefix, seq_along(valid_feature_group))
  names(valid_feature_id_map) <- sort(valid_feature_group)
  # invalid nuc
  invalid_feature_group <- unique(roi_node_df$feature_group) %>% subset(., !(. %in% valid_feature_group))
  invalid_feature_id_map <- paste0(invalid_feature_prefix, seq_along(invalid_feature_group))
  names(invalid_feature_id_map) <- invalid_feature_group
  
  ## Assign feature_id
  feature_id_map <- c(valid_feature_id_map, invalid_feature_id_map)
  roi_node_df <- mutate(roi_node_df, feature_id = feature_id_map[as.character(feature_group)])
  
  # # Visualizing graph
  # g <- tbl_graph(nodes=roi_node_df, edges=roi_edge_df, directed=FALSE)
  # g %>%
  #   activate(nodes) %>%
  #   mutate(feature_id = factor(feature_id)) %>%
  #   ggraph(layout = 'kk') +
  #   geom_edge_link() +
  #   geom_node_point(aes(color=feature_id), size = 1)
  
  # Preparing output ------------------------------------------------------------------------------
  out_df <- left_join(roi_pg_df, dplyr::select(roi_node_df, roi=roi_id, feature_id), by="roi", suffix=c("", "_x"))
  if(nrow(roi_fail_df) > 0){
    # Merged back with the held out ROIs
    out_df <- add_row(out_df, roi_fail_df)
  }
  return(out_df)
}


