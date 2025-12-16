# srun --job-name "interactive" --nodes=1 --cpus-per-task 16 --mem 128GB --time 12:00:00 --pty bash
# mamba activate R4.3_diffbind_chipseeker

# General packages
# library(yaml)
library(plyr)
library(dplyr)
library(stringr)
library(glue)
library(tidyr)
library(readxl)
library(writexl)

# genomics
library(DiffBind)
library(rtracklayer)
library(BiocIO)
library(csaw)
library(Rsamtools)
library(DESeq2)

# Ploting
library(ggrepel)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(circlize)
library(viridis)
library(RColorBrewer)
library(ComplexHeatmap)


#region Setup
proj_dir <- "." # Assume to be the root of this GitHub repo
setwd(proj_dir)
# source(paste0(proj_dir, "/R_session_setup.r"))
# .libPaths()

source(paste0(proj_dir, "/scripts/count_mat_2_pca.r"))

wd <- paste0(proj_dir, "/run_diffbind")
out_dir <- paste0(wd, "/data")
plot_dir <- paste0(wd, "/plot")
dir.create(out_dir, showWarnings=FALSE)
dir.create(plot_dir, showWarnings=FALSE)

mm10_blacklist_path <- "/path/to/ref/genome/mm10/blacklist_mm10_ENCFF547MET_noSelfOverlap_wChrM.bed"
mm10_blacklist <- BiocIO::import(mm10_blacklist_path)

# Read input --------------------------------------------------------------------------------------
force_run=FALSE

run_prefix <- "H3K27ac_NKKD" # Picking which analysis we want to run

if(run_prefix == "H3K27ac_NKKD"){
  samples_df <- read.csv(paste0(wd, "/input_tables/diffbind_H3K27ac_NKKD.csv"))
  ctrl_condition <- "NT_ctrl"
  
}else if(run_prefix == "ATAC_NKKD"){
  samples_df <- read.csv(paste0(wd, "/input_tables/diffbind_ATAC_NKKD.csv"))
  ctrl_condition <- "NT_ctrl"
  
}else if(run_prefix == NR5A2_Klf5KD_T7dedup){
  samples_df <- read.csv(paste0(wd, "/input_tables/diffbind_nr5a2_klf5KD_T7dedup.csv"))
  ctrl_condition = "NT_Ctrl"
  
}else{
  stop("Unknown run prefix ", run_prefix)
}


# Make sure the path is absolute
path_fields = c("bamReads", "bamControl", "Peaks")
for(i in seq_along(path_fields)){
  x <- samples_df[[path_fields[i]]]
  for(f in seq_along(x)){
    x[f] <- ifelse((is.na(x[f]) || x[f]==""), yes=NA, no=tools::file_path_as_absolute(x[f]))
  }
  samples_df[[path_fields[i]]] <- x
}


#region DBA count
# H3K27ac -----------------------------------------------------------------------------------------
dba_count_path <- paste0(out_dir, "/", run_prefix, "_dba-count.rds")
if(!file.exists(dba_count_path) || force_run){
  dba <- dba(sampleSheet=samples_df)
  dba <- dba.blacklist(dba, blacklist=mm10_blacklist)
  plot(dba)
  dba <- dba.count(dba) # Could take a while
  
  ## Normalization
  # dba <- dba.normalize(dba)
  dba <- dba.normalize(dba, normalize=DBA_NORM_LIB, 
                       method=DBA_DESEQ2, 
                       library=DBA_LIBSIZE_FULL, 
                       spikein="Spikein" %in% colnames(dba$samples))
  dba_norm_info <- dba.normalize(dba, bRetrieve=TRUE) # Retrieve normalize factor
  
  # dba.save(dba, file=str_remove(tools::file_path_sans_ext(basename(dba_count_path)), "^dba_"), dir=dirname(dba_count_path), ext=tools::file_ext(dba_count_path))
  saveRDS(dba, file=dba_count_path) # The read counting can take a while, back up the data
  
}else{
  cat("Loading existing DBA file (counts) ...")
  # dba <- dba.load(file=str_remove(tools::file_path_sans_ext(basename(dba_count_path)), "^dba_"), dir=dirname(dba_count_path), ext=tools::file_ext(dba_count_path))
  dba <- readRDS(dba_count_path)
  cat("Done\n")
}


#region DBA Diff
# Perform differential binding/accessibility
# this only works if we have >3 samples per condition
dba_diff_path <- paste0(out_dir, "/", run_prefix, "_dba-diff.rds")
if(!file.exists(dba_diff_path)){
  # dba <- dba.contrast(dba, design="~Condition", reorderMeta=list(Condition="NT_ctrl"))
  dba <- dba.contrast(dba, reorderMeta = list(Condition = ctrl_condition))
  dba <- dba.analyze(dba)
  
  saveRDS(dba, file=dba_diff_path)
  
}else{
  cat("Loading existing DBA file (diff binding) ...")
  dba <- readRDS(dba_diff_path)
  cat("Done\n")
}

# Extract information
dba_info <- dba.show(dba)
dba_info$peak_reads <- round(dba_info$Reads * dba_info$FRiP)
n_contrast <- length(dba$contrasts)

# Export merged peaks
merged_peaks <- as.data.frame(dba$merged)
colnames(merged_peaks) <- tolower(colnames(merged_peaks))
merged_peaks$chr <- dba$chrmap[dba$merged[ , 1]]
merged_peaks_gr <- makeGRangesFromDataFrame(merged_peaks)
BiocIO::export(merged_peaks_gr, con=paste0(out_dir, "/", run_prefix, "_DiffBind_merged_peaks.bed"))


#region Manual DESeq2
source(paste0(proj_dir, "/scripts/diffBind/dba_get_merged_peak.r"))
source(paste0(proj_dir, "/scripts/diffBind/dba_to_dds.r"))
source(paste0(proj_dir, "/scripts/diffBind/call_DBA.r"))

# MARK: Run DESeq2
manual_dds_path <- paste0(out_dir, "/manual_Diffbind_", run_prefix, "_dds.rds")
if(!file.exists(manual_dds_path)){
  dds <- dba_to_dds(dba)
  dds <- DESeq(dds)
  saveRDS(dds, manual_dds_path)
}else{
  cat("Reading dds ... ")
  dds <- readRDS(manual_dds_path)
  cat("Done\n")
}


# # doing quick check
# head(DESeq2::counts(dds, normalized=TRUE))
# head(dba$binding)
# colData(dds)$Condition

## Set up contrast list
if(run_prefix %in% "NR5A2_Klf5KD_T7dedup"){
  # Nr5a2 binding in Klf5 KD
  contrast_list <- list(
    Klf5_KD_vs_NT_Ctrl = c("Condition", "Klf5_KD", "NT_Ctrl")
  )
}else if(run_prefix %in% c("H3K27ac_NKKD", "ATAC_NKKD")){
  # For Atac and H3K27ac
  contrast_list <- list(
    Nr5a2_KD_vs_NT_ctrl = c("Condition", "Nr5a2_KD", "NT_ctrl"),
    Klf5_KD_vs_NT_ctrl = c("Condition", "Klf5_KD", "NT_ctrl"),
    dKD_vs_NT_ctrl = c("Condition", "dKD", "NT_ctrl")
  )
}else{
  stop("Unknown run_prefix: ", run_prefix)
}

# Set contrast name
for(i in seq_along(contrast_list)){
  names(contrast_list)[i] <- paste0(contrast_list[[i]][2], "_vs_", contrast_list[[i]][3])
}


# MARK: plot (DBA diff)
DBA_colmap <- c(Gain = "orange", No_Change = "gray", Lost = "darkorchid4")
tmp_plot_DBA <- function(res_df, plot_title="Differential binding analysis (DBA)", subtitle="", DBA_colmap=c(Gain = "orange", No_Change = "gray", Lost = "darkorchid4")){
  p <- res_df %>%
    dplyr::mutate(DBA_factor = factor(DBA, level = c("No_Change", "Lost", "Gain"))) %>%
    dplyr::arrange(DBA_factor) %>%
    ggplot(aes(x = baseMean, y = log2FoldChange, color = DBA)) +
    geom_hline(yintercept = 0) +
    geom_point(alpha = 0.8) +
    scale_color_manual(values = DBA_colmap) +
    labs(title = plot_title, subtitle = subtitle) +
    scale_x_log10() +
    theme_classic() +
    theme(legend.position="none")
  return(p)
}

tmp_resDf_2_bed <- function(res_df, DBA_colmap = c(Gain = "orange", No_Change = "gray", Lost = "darkorchid4")){
  DBA_rgb_mat <- col2rgb(DBA_colmap[res_df$DBA])
  res_df$rgb <- apply(DBA_rgb_mat, 2, paste0, collapse = ",")
  all_region_bed <- res_df %>%
    tibble::rownames_to_column(var = "id") %>%
    dplyr::mutate(
      name = paste0(id, "_", DBA),
      strand = "*",
      t_start = start,
      t_end = end
    ) %>%
    dplyr::select(seqnames, start, end, name, score = log2FoldChange, strand, t_start, t_end, rgb, DBA) %>%
    dplyr::arrange(seqnames, start, end)
  
  return(all_region_bed)
}

plot_suffix <- ""; if("MA_ylim" %in% ls()){rm(MA_ylim)}
MA_ylim = c(-7, 7); plot_suffix = "_fixed_ylim"
# RColorBrewer::display.brewer.all()
res_list <- list()
for(i in seq_along(contrast_list)){
  # For DBA diff
  contrast_name <- names(contrast_list)[i]
  contrast_vec <- contrast_list[[i]]

  message("Checking results from ", contrast_name)
  #TODO: Adding region location to the res_df
  region_cols <- c("seqnames", "start", "end")
  if ("DESeq2" %in% names(dba)) {
    message("  Ploting result from DBA")
    dba_contrast_vec <- contrast_vec
    dba_contrast_vec[1] <- "Treatment" # It use this column of colData for some reason
    res_df <- as_tibble(DESeq2::results(dba$DESeq2$DEdata, contrast = dba_contrast_vec, alpha = 0.05))
    res_df <- call_DBA(res_df)
    deg_stats_txt <- as.data.frame(table(res_df$DBA)) %>%
      dplyr::mutate(txt = paste0(Var1, ": ", Freq)) %>%
      dplyr::pull(txt) %>%
      paste0(collapse = " | ")
    cur_stt <- paste0(
      "Comparison: ", dba_contrast_vec[2], " vs ", dba_contrast_vec[3],
      "\nNorm method: ", dba$norm$DESeq2$norm.calc,
      "\nLibrary calculation method: ", dba$norm$DESeq2$lib.calc,
      "\nDBA: ", deg_stats_txt
    )
    ylim_vec <- c(floor(min(res_df$log2FoldChange)), ceiling(max(res_df$log2FoldChange)))
    y_break_max = max(abs(res_df$log2FoldChange))
    if (("MA_ylim" %in% ls())) {
      ylim_vec <- MA_ylim
      y_break_max <- max(c(y_break_max, abs(MA_ylim)))
    }
    if (y_break_max %% 2 != 0) y_break_max <- y_break_max + 1
    y_break_vec <- seq(-y_break_max, y_break_max, by = 2)
    p1 <- tmp_plot_DBA(res_df, plot_title = paste0("DBA | ", run_prefix), subtitle = cur_stt)
    p1 <- p1 + scale_y_continuous(breaks = y_break_vec, limits = ylim_vec)
    ggsave(
      plot = p1, filename = paste0(plot_dir, "/DiffBind_MA_", run_prefix, "_", names(contrast_list)[i], plot_suffix, ".pdf"),
      dpi = 300, units = "in", width = 5, height = 5
    )
    
    region_df <- as.data.frame(dba$merged)
    colnames(region_df) <- c("seqnames", "start", "end")
    region_df$seqnames <- dba$chrmap[region_df$seqnames]
    res_df <- res_df %>%
      dplyr::select(-dplyr::any_of(region_cols)) %>%
      dplyr::bind_cols(region_df) %>%
      dplyr::select(dplyr::all_of(region_cols), dplyr::everything())
    
    all_region_bed <- tmp_resDf_2_bed(res_df)
    write.table(dplyr::select(all_region_bed, -DBA),
      file = paste0(out_dir, "/DiffBind_", run_prefix, "_", names(contrast_list)[i], "_all_regions.bed"),
      sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
    )
    all_region_bed %>%
      dplyr::filter(str_detect(tolower(DBA), "(gain)|(los[st])")) %>%
      write.table(
        file = paste0(out_dir, "/DiffBind_", run_prefix, "_", names(contrast_list)[i], "_DBA_only.bed"),
        sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
      )
    res_list[[paste0(names(contrast_list)[i], "_default")]] <- res_df
  }
  
  # For manual DESeq2
  message("  Ploting result from manual DESeq2")
  res_df <- results(dds, contrast = contrast_list[[i]], alpha = 0.05, lfcThreshold = 0)
  res_df <- call_DBA(as.data.frame(res_df))
  # Add regions to res_df
  if(any(!region_cols %in% colnames(res_df))){
    region_df <- as.data.frame((rowRanges(dds)))
    region_df <- dplyr::select(region_df, dplyr::all_of(region_cols))
    res_df <- res_df %>%
      dplyr::select(-dplyr::any_of(region_cols)) %>% 
      dplyr::bind_cols(region_df[rownames(res_df), , drop = FALSE]) %>%
      dplyr::select(dplyr::all_of(region_cols), dplyr::everything())
  }
  deg_stats_txt <- as.data.frame(table(res_df$DBA)) %>%
    dplyr::mutate(txt = paste0(Var1, ": ", Freq)) %>%
    dplyr::pull(txt) %>%
    paste0(collapse = " | ")
  cur_stt <- paste0(
    "Comparison: ", contrast_vec[2], " vs ", contrast_vec[3],
    "\nNorm method: ", dba$norm$DESeq2$norm.calc,
    "\nLibrary calculation method: ", dba$norm$DESeq2$lib.calc,
    "\nDBA: ", deg_stats_txt
  )
  ylim_vec <- c(floor(min(res_df$log2FoldChange)), ceiling(max(res_df$log2FoldChange)))
  y_break_max = max(abs(res_df$log2FoldChange))
  if (("MA_ylim" %in% ls())) {
    ylim_vec <- MA_ylim
    y_break_max <- max(c(y_break_max, abs(MA_ylim)))
  }
  if (y_break_max %% 2 != 0) y_break_max <- y_break_max + 1
  y_break_vec <- seq(-y_break_max, y_break_max, by = 2)
  p <- tmp_plot_DBA(res_df, plot_title = paste0("DBA | (manual) ", run_prefix), subtitle = cur_stt)
  p <- p + scale_y_continuous(breaks = y_break_vec, limits = ylim_vec)
  ggsave(
    plot = p, filename = paste0(plot_dir, "/manual_DiffBind_MA_", run_prefix, "_", names(contrast_list)[i], plot_suffix, ".pdf"),
    dpi = 300, units = "in", width = 5, height = 5
  )
  # Export DBA regions
  all_region_bed <- tmp_resDf_2_bed(res_df)
  write.table(dplyr::select(all_region_bed, -DBA),
    file = paste0(out_dir, "/manual_DiffBind_", run_prefix, "_", names(contrast_list)[i], "_all_regions.bed"),
    sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
  )

  all_region_bed %>%
    dplyr::filter(str_detect(tolower(DBA), "(gain)|(los[st])")) %>%
    write.table(
      file = paste0(out_dir, "/manual_DiffBind_", run_prefix, "_", names(contrast_list)[i], "_DBA_only.bed"),
      sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
    )

  # Concat to the main list
  res_list[[names(contrast_list)[i]]] <- res_df
}


## Save res_list ----------------------------------------------------------------------------------
saveRDS(res_list, file=paste0(out_dir, "/DiffBind_", run_prefix, "_results_list.rds"))
writexl::write_xlsx(res_list, path=paste0(out_dir, "/DiffBind_", run_prefix, "_results.xlsx"))

# Results list for DBA only
res_list_dba <- list()
for(i in seq_along(res_list)){
  df <- res_list[[i]] %>% 
    dplyr::filter(DBA != "No_Change")
  if(nrow(df) > 0){
    res_list_dba[[names(res_list)[i]]] <- df
  }
}

if(length(res_list_dba) > 0){
  write_xlsx(res_list_dba, path=paste0(out_dir, "/DiffBind_", run_prefix, "_results_DBA_only.xlsx"))
}



cor_colmap <- circlize::colorRamp2(
  breaks = seq(0.5, 1, length.out = 8),
  colors = rev(RColorBrewer::brewer.pal(name = "RdGy", n = 8))
)


# Plot correlation
pcc_mat <- cor(DESeq2::counts(dds, normalized = TRUE), method = "pearson")
scc_mat <- cor(DESeq2::counts(dds, normalized = TRUE), method = "spearman")
hm_pcc <- ComplexHeatmap::Heatmap(pcc_mat, col = cor_colmap, name = "PCC")
hm_scc <- ComplexHeatmap::Heatmap(scc_mat, col = cor_colmap, name = "SCC")
# Save heatmap plot
hm_basename <- paste0("manual_DiffBind_Correlation_", run_prefix)
{
  pdf(paste0(plot_dir, "/", hm_basename, ".pdf"), width = 8, height = 6)
  ComplexHeatmap::draw(hm_pcc)
  ComplexHeatmap::draw(hm_scc)
  dev.off()
}




