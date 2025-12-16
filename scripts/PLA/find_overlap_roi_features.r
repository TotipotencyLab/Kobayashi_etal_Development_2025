#' @description
#' This function taks two previously defined cellular features and determined whether they are overlap with each other or not.
#' 
#' 
#' @param feature_df_1 data.frame of cellular feature (e.g., the output of `define_feature_group` function)
#' @param feature_df_2 another cellular feature dataframe
#' @param feature_1_regex Character regular expression for filtering (by feature name) which features from feature_df_1 to be uased
#' @param feature_2_regex Character Same as `feature_2_regex`
#' @param min_intersect_ratio Numeric. Value between 0 to 1. Minimum ratio of ROI intersection to be considered as overlap
#' @param min_ratio_roi_overlap Numeric. Value between 0 to 1. Minimum ratio of ROI that are overlap with each other to consider the two features as overlap.
#'  ratio_roi_overlap = length(overlap_roi) / min(length(ROI_1), length(ROI_2))
#'  if 1, all ROI in the feature that has smaller number of ROI has to be part of the overlap

# min_share_z_ratio:
#   share_z_ratio = length(share_z with overlap) / min(length(z_1), length(z_2))
#   if 1, feature with smaller number of z-stack spaning must be fully overlap with the larger z-stack spanning feature.


find_overlap_roi_features <- function(feature_df_1, feature_df_2, 
                                      feature_1_regex="feature_", feature_2_regex="feature_",
                                      min_intersect_ratio=0.5, 
                                      # min_share_z_ratio = 1,
                                      min_ratio_roi_overlap=1){
  
  # Detect feature 
  F1_df <- dplyr::filter(feature_df_1, str_detect(feature_id, feature_1_regex))
  F2_df <- dplyr::filter(feature_df_2, str_detect(feature_id, feature_2_regex))
  
  if(nrow(F1_df)==0){stop("Couldn't find feature 1")}
  if(nrow(F2_df)==0){stop("Couldn't find feature 2")}
  
  # F1_feature <- unique(F1_df$feature_id)
  # F2_feature <- unique(F2_df$feature_id)
  # roi_feature_map_1 <- dplyr::pull(F1_df, feature_id, name=roi)
  # roi_feature_map_2 <- dplyr::pull(F2_df, feature_id, name=roi)
  
  # # Area
  F1_df$area <- st_area(F1_df$geometry)
  F2_df$area <- st_area(F2_df$geometry)
  # roi_area_map_1 <- dplyr::pull(F1_df, area, name=roi)
  # roi_area_map_2 <- dplyr::pull(F2_df, area, name=roi)
  
  # # Comparing features
  # z_range_list_1 <- split(F1_df$z, F1_df$feature_id)
  # z_range_list_2 <- split(F2_df$z, F2_df$feature_id)
  # roi_z_map_1 <- dplyr::pull(F1_df, z, name=roi)
  # roi_z_map_2 <- dplyr::pull(F2_df, z, name=roi)
  # 
  # roi_list_1 <- split(F1_df$roi, F1_df$feature_id)
  # roi_list_2 <- split(F2_df$roi, F2_df$feature_id)
  
  # Find overlap of all ROIs ----------------------------------------------------------------------
  roi_ovl_mat <- st_intersects(F1_df$geometry, F2_df$geometry, sparse=F)
  rownames(roi_ovl_mat) <- F1_df$roi
  colnames(roi_ovl_mat) <- F2_df$roi
  
  roi_ovl_pair_df <- as.data.frame(roi_ovl_mat) %>% 
    rownames_to_column(var="roi_1") %>% 
    tidyr::gather(key="roi_2", value="overlap", -roi_1) %>% 
    as_tibble()
  
  roi_ovl_pair_df <- roi_ovl_pair_df %>% 
    left_join(., dplyr::select(F1_df, roi_1=roi, z_1=z, area_1=area, f_id_1=feature_id), by="roi_1") %>% 
    left_join(., dplyr::select(F2_df, roi_2=roi, z_2=z, area_2=area, f_id_2=feature_id), by="roi_2") %>% 
    # Rearrange (for debugging)
    dplyr::select(starts_with("roi"), overlap, starts_with("f_id"), starts_with("z"), starts_with("area"), everything()) %>% 
    mutate(min_area = pmin(area_1, area_2),
           int_area = 0,
           ovl_idx = paste0("idx", seq_along(roi_1)))
  
  # Only consider roi on the same z-stack
  roi_ovl_pair_df <- roi_ovl_pair_df %>% 
    mutate(overlap = (overlap & (z_1 == z_2)))
  
  # Calculate overlap area ------------------------------------------------------------------------
  # Only do this for those with known to be overlap
  if(sum(roi_ovl_pair_df$overlap)>0){
    roi_ovl_pair_df_fil <- dplyr::filter(roi_ovl_pair_df, overlap)
    for(i in 1:nrow(roi_ovl_pair_df_fil)){
      pg1 <- F1_df$geometry[F1_df$roi == roi_ovl_pair_df_fil$roi_1[i]]
      pg2 <- F2_df$geometry[F2_df$roi == roi_ovl_pair_df_fil$roi_2[i]]
      # plot(c(pg1, pg2)) # For visual inspection
      ovl_area <- st_area(st_intersection(pg1, pg2))
      
      roi_ovl_pair_df_fil$int_area[i] <- ovl_area
    }
    
    # Merge the data back
    roi_ovl_pair_df <- left_join(dplyr::select(roi_ovl_pair_df, -int_area), 
                                 dplyr::select(roi_ovl_pair_df_fil, ovl_idx, int_area), by="ovl_idx") %>% 
      replace_na(replace=list(int_area=0))
    
    # Calculate intersect ratio
    roi_ovl_pair_df <- roi_ovl_pair_df %>% 
      mutate(int_ratio = as.double(int_area/min_area),
             overlap = (overlap & (int_ratio >= min_intersect_ratio)))
    
    # # For checking
    # roi_ovl_pair_df %>% 
    #   dplyr::filter(int_ratio != 0) %>% 
    #   dplyr::pull(int_ratio) %>% 
    #   range()
    
    # Detecting overlap features --------------------------------------------------------------------
    # Summarizing overlap by feature
    feature_pair_stats <- roi_ovl_pair_df %>% 
      group_by(f_id_1, f_id_2) %>% 
      reframe(
        # # Z-filtering
        # n_z_1 = length(unique(z_1)),
        # n_z_2 = length(unique(z_2)),
        # min_z = pmin(n_z_1, n_z_2),
        # n_z_1_ovl = length(unique(z_1[overlap])),
        # n_z_2_ovl = length(unique(z_2[overlap])),
        # ratio_z_ovl = pmin(n_z_1_ovl, n_z_2_ovl) / min_z,
        
        ## ROI ovl ratio filtering
        n_roi_1 = length(unique(roi_1)),
        n_roi_2 = length(unique(roi_2)),
        min_roi = pmin(n_roi_1, n_roi_2),
        n_roi_1_ovl = length(unique(roi_1[overlap])),
        n_roi_2_ovl = length(unique(roi_2[overlap])),
        ratio_roi_ovl = pmin(n_roi_1_ovl, n_roi_2_ovl) / min_roi
      )
    
    # Preparing output ------------------------------------------------------------------------------
    feature_ovl_df <- feature_pair_stats %>% 
      dplyr::filter(ratio_roi_ovl >= min_ratio_roi_overlap) %>% 
      dplyr::select(feature_id_1 = f_id_1,
                    feature_id_2 = f_id_2)
    
  }else{
    # i.e., no overlap to begin with
    feature_ovl_df <- tibble()
  }
  
  return(feature_ovl_df)
}


