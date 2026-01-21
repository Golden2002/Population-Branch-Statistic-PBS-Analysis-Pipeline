#!/usr/bin/env Rscript
# ============================================================================
# PBS (Population Branch Statistic) Analysis Script
# Function: Calculate PBS branch statistics based on pairwise Fst data to identify selection signals in target populations
# Author: Developed based on user requirements
# Date: 2024
# ============================================================================

library(dplyr)
library(data.table)
library(tidyr)
library(purrr)

# ============================================================================
# 1. Configuration and Parameter Settings
# ============================================================================

# Set working directory (modify according to actual situation)
setwd("PATH_TO_YOUR_WORKING_DIRECTORY")

# Output directory
output_dir <- "PBS_analysis_results"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Define experimental design list: each design is a triple (target, ref1, ref2)
# Additional designs can be added as needed
designs <- list(
  "Target_Ref1_Ref2" = c("Target", "Reference1", "Reference2"),
  "Target_Ref3_Ref4" = c("Target", "Reference3", "Reference4")
  # Example: add more designs
  # "Ref1_Ref3_Ref4" = c("Reference1", "Reference3", "Reference4")
)

# Define file naming pattern (adjust according to actual filenames)
# Assumed filename format: Target_vs_Reference.windowed.weir.fst
file_pattern <- "_vs_"

# PBS calculation parameters
top_percent <- 0.01  # Take top 1%
min_fst <- 1e-6      # Minimum Fst value to avoid log calculation issues

# ============================================================================
# 2. Data Reading Functions
# ============================================================================

#' Extract population pair information from filename
#' 
#' @param filename Filename
#' @param pattern Separator pattern
#' @return Vector containing two population names
extract_pop_pair <- function(filename, pattern = "_vs_") {
  # Remove file extension
  base_name <- tools::file_path_sans_ext(basename(filename))
  base_name <- sub("\\.windowed\\.weir\\.fst$", "", base_name)
  
  # Split population names
  pops <- unlist(strsplit(base_name, pattern))
  if (length(pops) != 2) {
    warning(paste("Filename", filename, "does not match expected population pair format"))
    return(NULL)
  }
  return(pops)
}

#' Read Fst file and standardize
#' 
#' @param file_path Fst file path
#' @return Data frame containing chromosome positions and Fst values
# Create a function to standardize pair names
standardize_pair <- function(pop1, pop2) {
  # Iterate through designs to find matches
  for(design in designs) {
    if(pop1 %in% design && pop2 %in% design) {
      # Sort according to order in design
      # Ensure returned order matches the order in design
      ordered <- design[design %in% c(pop1, pop2)]
      return(paste(ordered, collapse = "_vs_"))
    }
  }
  # No match found, sort alphabetically
  return(paste(sort(c(pop1, pop2)), collapse = "_vs_"))
}

read_fst_file <- function(file_path) {
  cat("Reading file:", basename(file_path), "\n")
  
  # Extract population pair information
  pops <- extract_pop_pair(file_path, file_pattern)

  # Remove residual .windowed.weir
  pops <- gsub("\\.windowed\\.weir$", "", pops)

  if (is.null(pops)) return(NULL)
  
  # Debug information
  cat("pops extracted:", pops, "\n")
  cat("pair generated:", standardize_pair(pops[1], pops[2]), "\n")
  
  # Read data
  fst_data <- tryCatch({
    fread(file_path,
          col.names = c("chr", "start", "end", "n_snp", "mean_fst", "weighted_fst"),
          select = 1:6)
  }, error = function(e) {
    warning(paste("Unable to read file:", file_path, "Error:", e$message))
    return(NULL)
  })
  
  if (is.null(fst_data) || nrow(fst_data) == 0) {
    warning(paste("File is empty or reading failed:", file_path))
    return(NULL)
  }
  
  # Data cleaning
  fst_data <- fst_data %>%
    filter(!is.na(weighted_fst)) %>%
    mutate(
      weighted_fst = ifelse(weighted_fst < 0, 0, weighted_fst),
      weighted_fst = ifelse(weighted_fst >= 1, 1 - 1e-6, weighted_fst),
      pop1 = pops[1],
      pop2 = pops[2],
      pair = gsub("\\.windowed\\.weir$", "", standardize_pair(pops[1], pops[2]))  # Remove suffix
    ) %>%
    ungroup() %>%
    select(chr, start, end, weighted_fst, pair, pop1, pop2)
  
  return(fst_data)
}

#' Load all Fst files
#' 
#' @param data_dir Data directory
#' @param pattern File pattern
#' @return List containing all Fst data
load_all_fst_data <- function(data_dir = ".", pattern = "\\.windowed\\.weir\\.fst$") {
  # Find all Fst files
  fst_files <- list.files(data_dir, pattern = pattern, full.names = TRUE, recursive = TRUE)
  
  if (length(fst_files) == 0) {
    stop("No Fst files found, please check directory and file pattern")
  }
  
  cat("Found", length(fst_files), "Fst files\n")
  
  # Read all files
  all_data <- list()
  for (file in fst_files) {
    data <- read_fst_file(file)
    if (!is.null(data)) {
      # Use standardized pair names as list keys
      all_data[[unique(data$pair)]] <- data
    }
  }
  
  cat("Successfully loaded", length(all_data), "valid Fst datasets\n")
  return(all_data)
}

# ============================================================================
# 3. PBS Calculation Functions
# ============================================================================

#' Calculate branch length T = -log(1 - Fst)
#' 
#' @param fst Fst value vector
#' @param min_fst Minimum Fst value to avoid numerical issues
#' @return Branch length T value vector
calculate_T <- function(fst, min_fst = 1e-6) {
  # Ensure Fst is within valid range
  fst_adj <- pmax(fst, min_fst)
  fst_adj <- pmin(fst_adj, 1 - min_fst)
  
  # Calculate T = -log(1 - Fst)
  T_values <- -log(1 - fst_adj)
  return(T_values)
}

#' Calculate PBS values
#' 
#' @param fst_data List containing all Fst data
#' @param target Target population
#' @param ref1 Reference population 1
#' @param ref2 Reference population 2
#' @param min_fst Minimum Fst value
#' @return Data frame containing PBS values
calculate_PBS <- function(fst_data, target, ref1, ref2, min_fst = 1e-6) {
  cat("Calculating PBS:", target, "|", ref1, "|", ref2, "\n")
  
  # Use standardize_pair function to construct keys for three population pairs
  pairs <- list(
    target_ref1 = standardize_pair(target, ref1),
    target_ref2 = standardize_pair(target, ref2),
    ref1_ref2 = standardize_pair(ref1, ref2)
  )
  
  # Debug information: print keys being searched
  cat("Searching for data with the following keys:\n")
  for (pair_name in names(pairs)) {
    cat("  ", pair_name, ":", pairs[[pair_name]], "\n")
  }
  
  # Check available keys
  cat("Available data keys:\n")
  print(names(fst_data))
  
  # Check if required data exists
  missing_pairs <- sapply(pairs, function(p) is.null(fst_data[[p]]))
  if (any(missing_pairs)) {
    missing_names <- names(pairs)[missing_pairs]
    warning(paste("Missing Fst data for the following population pairs:", paste(missing_names, collapse = ", ")))
    return(NULL)
  }
  
  # Extract Fst data for three population pairs
  df_target_ref1 <- fst_data[[pairs$target_ref1]] %>%
    select(chr, start, end, weighted_fst) %>%
    rename(fst_TR1 = weighted_fst)
  
  df_target_ref2 <- fst_data[[pairs$target_ref2]] %>%
    select(chr, start, end, weighted_fst) %>%
    rename(fst_TR2 = weighted_fst)
  
  df_ref1_ref2 <- fst_data[[pairs$ref1_ref2]] %>%
    select(chr, start, end, weighted_fst) %>%
    rename(fst_R1R2 = weighted_fst)
  
  # Merge three data frames
  combined <- df_target_ref1 %>%
    inner_join(df_target_ref2, by = c("chr", "start", "end")) %>%
    inner_join(df_ref1_ref2, by = c("chr", "start", "end"))
  
  if (nrow(combined) == 0) {
    warning("Three population pairs have no common windows")
    return(NULL)
  }
  
  # Calculate branch length T
  combined <- combined %>%
    mutate(
      T_TR1 = calculate_T(fst_TR1, min_fst),
      T_TR2 = calculate_T(fst_TR2, min_fst),
      T_R1R2 = calculate_T(fst_R1R2, min_fst)
    )
  
  # Calculate PBS: PBS = (T_TR1 + T_TR2 - T_R1R2) / 2
  combined <- combined %>%
    mutate(
      PBS = (T_TR1 + T_TR2 - T_R1R2) / 2,
      # Ensure PBS is non-negative
      PBS = pmax(PBS, 0),
      design = paste(target, ref1, ref2, sep = "_")
    ) %>%
    select(chr, start, end, fst_TR1, fst_TR2, fst_R1R2, T_TR1, T_TR2, T_R1R2, PBS, design)
  
  cat("Successfully calculated PBS values for", nrow(combined), "windows\n")
  return(combined)
}

# ============================================================================
# 4. Analysis and Result Processing Functions
# ============================================================================

#' Extract top percentage results
#' 
#' @param pbs_data PBS data frame
#' @param percent Top percentage (between 0-1)
#' @return Data frame of top percentage
extract_top_percent <- function(pbs_data, percent = 0.01) {
  if (percent <= 0 || percent >= 1) {
    stop("Percentage must be between 0 and 1")
  }
  
  # Calculate threshold
  n_sites <- nrow(pbs_data)
  n_top <- ceiling(n_sites * percent)
  
  # Sort by PBS value in descending order, take top percentage
  top_sites <- pbs_data %>%
    arrange(desc(PBS)) %>%
    head(n_top)
  
  cat("Extracting top", percent*100, "% (", n_top, "sites)\n", sep = "")
  return(top_sites)
}

#' Calculate intersection of multiple designs
#' 
#' @param top_results_list List containing top results of multiple designs
#' @return Data frame of intersection sites
calculate_intersection <- function(top_results_list) {
  if (length(top_results_list) < 2) {
    warning("At least two designs are required to calculate intersection")
    return(NULL)
  }
  
  # Extract site identifiers for each design
  site_sets <- lapply(names(top_results_list), function(design_name) {
    df <- top_results_list[[design_name]]
    df %>%
      mutate(site_id = paste(chr, start, end, sep = ":")) %>%
      pull(site_id)
  })
  names(site_sets) <- names(top_results_list)
  
  # Calculate intersection
  intersect_sites <- Reduce(intersect, site_sets)
  
  cat("Number of intersection sites:", length(intersect_sites), "\n")
  
  # Create intersection data frame
  if (length(intersect_sites) > 0) {
    # Extract complete information of intersection sites from first design
    first_design <- names(top_results_list)[1]
    first_df <- top_results_list[[first_design]] %>%
      mutate(site_id = paste(chr, start, end, sep = ":")) %>%
      filter(site_id %in% intersect_sites) %>%
      select(-site_id)
    
    # Add PBS values from all designs
    intersect_df <- first_df
    for (design_name in names(top_results_list)[-1]) {
      design_df <- top_results_list[[design_name]] %>%
        mutate(site_id = paste(chr, start, end, sep = ":")) %>%
        filter(site_id %in% intersect_sites) %>%
        select(site_id, PBS)
      
      col_name <- paste0("PBS_", design_name)
      intersect_df <- intersect_df %>%
        mutate(site_id = paste(chr, start, end, sep = ":")) %>%
        left_join(design_df %>% rename(!!col_name := PBS), by = "site_id") %>%
        select(-site_id)
    }
    
    return(intersect_df)
  } else {
    return(NULL)
  }
}

#' Save results
#' 
#' @param data Data frame to save
#' @param filename Filename
#' @param output_dir Output directory
#' @param sort_by Sorting method: "PBS" or "position"
save_results <- function(data, filename, output_dir, sort_by = "PBS") {
  if (is.null(data) || nrow(data) == 0) {
    cat("No data to save:", filename, "\n")
    return(FALSE)
  }
  
  # Sort according to sorting method
  if (sort_by == "PBS") {
    # Find columns containing "PBS" (may have multiple)
    pbs_cols <- grep("^PBS", colnames(data), value = TRUE)
    if (length(pbs_cols) > 0) {
      # Sort in descending order by first PBS column
      sorted_data <- data %>% arrange(desc(.[[pbs_cols[1]]]))
    } else {
      sorted_data <- data
    }
  } else if (sort_by == "position") {
    sorted_data <- data %>% arrange(chr, start, end)
  } else {
    sorted_data <- data
  }
  
  # Create complete file path
  filepath <- file.path(output_dir, filename)
  
  # Save as CSV
  fwrite(sorted_data, filepath, row.names = FALSE)
  cat("Results saved to:", filepath, "\n")
  
  # Also save as BED format (convenient for genome browser viewing)
  if ("chr" %in% colnames(sorted_data) && "start" %in% colnames(sorted_data) && "end" %in% colnames(sorted_data)) {
    bed_data <- sorted_data %>%
      select(chr, start, end) %>%
      mutate(chr = as.character(chr),
             start_bed = start - 1,  # BED format is 0-based
             name = paste("site", 1:n(), sep = "_")) %>%
      select(chr, start_bed, end, name)
    
    bed_filepath <- sub("\\.csv$", ".bed", filepath)
    fwrite(bed_data, bed_filepath, col.names = FALSE, sep = "\t", quote = FALSE)
    cat("BED format saved to:", bed_filepath, "\n")
  }
  
  return(TRUE)
}

# ============================================================================
# 5. Main Functions
# ============================================================================

#' Main analysis function
#' 
#' @param designs Experimental design list
#' @param top_percent Top percentage to extract
#' @param output_dir Output directory
run_PBS_analysis <- function(designs, top_percent = 0.01, output_dir = "PBS_results") {
  # Create output directory
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  cat("Starting PBS analysis...\n")
  cat("Number of experimental designs:", length(designs), "\n")
  
  # 1. Load all Fst data
  cat("\nStep 1: Loading Fst data...\n")
  fst_data <- load_all_fst_data()
  
  # 2. Calculate PBS for each design
  cat("\nStep 2: Calculating PBS for each design...\n")
  pbs_results <- list()
  
  for (design_name in names(designs)) {
    pops <- designs[[design_name]]
    if (length(pops) != 3) {
      warning(paste("Design", design_name, "must contain 3 populations"))
      next
    }
    
    pbs_df <- calculate_PBS(fst_data, pops[1], pops[2], pops[3], min_fst = min_fst)
    
    if (!is.null(pbs_df)) {
      pbs_results[[design_name]] <- pbs_df
      
      # Save complete PBS results
      full_filename <- paste0("PBS_full_", design_name, ".csv")
      save_results(pbs_df, full_filename, output_dir, sort_by = "PBS")
    }
  }
  
  if (length(pbs_results) == 0) {
    stop("No PBS values successfully calculated for any design")
  }
  
  # 3. Extract top percentage results for each design
  cat("\nStep 3: Extracting top", top_percent*100, "% results for each design...\n", sep = "")
  top_results <- list()
  
  for (design_name in names(pbs_results)) {
    top_df <- extract_top_percent(pbs_results[[design_name]], top_percent)
    top_results[[design_name]] <- top_df
    
    # Save top results (two sorting methods)
    top_filename_pbs <- paste0("PBS_top_", top_percent*100, "percent_", design_name, "_byPBS.csv")
    top_filename_pos <- paste0("PBS_top_", top_percent*100, "percent_", design_name, "_byPosition.csv")
    
    save_results(top_df, top_filename_pbs, output_dir, sort_by = "PBS")
    save_results(top_df, top_filename_pos, output_dir, sort_by = "position")
  }
  
  # 4. Calculate intersection
  cat("\nStep 4: Calculating intersection of top results across designs...\n")
  if (length(top_results) >= 2) {
    intersect_df <- calculate_intersection(top_results)
    
    if (!is.null(intersect_df) && nrow(intersect_df) > 0) {
      # Save intersection results (two sorting methods)
      intersect_filename_pbs <- paste0("PBS_intersection_top", top_percent*100, "percent_byPBS.csv")
      intersect_filename_pos <- paste0("PBS_intersection_top", top_percent*100, "percent_byPosition.csv")
      
      save_results(intersect_df, intersect_filename_pbs, output_dir, sort_by = "PBS")
      save_results(intersect_df, intersect_filename_pos, output_dir, sort_by = "position")
    } else {
      cat("No intersection sites found\n")
    }
  }
  
  # 5. Generate summary report
  cat("\nStep 5: Generating summary report...\n")
  generate_summary_report(pbs_results, top_results, output_dir, top_percent)
  
  cat("\nAnalysis completed! Results saved in:", output_dir, "\n")
  
  # Return result list
  return(list(
    pbs_results = pbs_results,
    top_results = top_results,
    intersect_results = if (exists("intersect_df")) intersect_df else NULL
  ))
}

#' Generate summary report
#' 
#' @param pbs_results Complete PBS results for each design
#' @param top_results Top results for each design
#' @param output_dir Output directory
#' @param top_percent Top percentage
generate_summary_report <- function(pbs_results, top_results, output_dir, top_percent) {
  summary_lines <- c(
    "PBS Analysis Summary Report",
    paste("Generation time:", Sys.time()),
    paste("Number of designs analyzed:", length(pbs_results)),
    paste("Top percentage:", top_percent*100, "%"),
    "\nStatistical information for each design:"
  )
  
  for (design_name in names(pbs_results)) {
    pbs_df <- pbs_results[[design_name]]
    top_df <- top_results[[design_name]]
    
    design_summary <- c(
      paste("\nDesign:", design_name),
      paste("  Total windows:", nrow(pbs_df)),
      paste("  Mean PBS:", round(mean(pbs_df$PBS, na.rm = TRUE), 6)),
      paste("  Median PBS:", round(median(pbs_df$PBS, na.rm = TRUE), 6)),
      paste("  Maximum PBS:", round(max(pbs_df$PBS, na.rm = TRUE), 6)),
      paste("  Top", top_percent*100, "% sites:", nrow(top_df)),
      paste("  Minimum PBS in top sites:", round(min(top_df$PBS, na.rm = TRUE), 6))
    )
    
    summary_lines <- c(summary_lines, design_summary)
  }
  
  # Write to file
  report_path <- file.path(output_dir, "PBS_analysis_summary.txt")
  writeLines(summary_lines, report_path)
  cat("Summary report saved to:", report_path, "\n")
}

# ============================================================================
# 6. Execute Analysis
# ============================================================================

cat(strrep("=", 60), "\n\n")
cat("PBS Branch Statistic Analysis Script\n")
cat(strrep("=", 60), "\n\n")

# Execute analysis
analysis_results <- tryCatch({
  run_PBS_analysis(designs, top_percent, output_dir)
}, error = function(e) {
  cat("Error during analysis:", e$message, "\n")
  NULL
})

# ============================================================================
# 7. Optional: Visualization Functions (if needed)
# ============================================================================

#' Plot PBS distribution
#' 
#' @param pbs_data PBS data frame
#' @param design_name Design name
#' @param output_dir Output directory
plot_PBS_distribution <- function(pbs_data, design_name, output_dir) {
  library(ggplot2)
  
  # Create plot directory
  plot_dir <- file.path(output_dir, "plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  # PBS distribution histogram
  p1 <- ggplot(pbs_data, aes(x = PBS)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = quantile(pbs_data$PBS, 0.99, na.rm = TRUE), 
               color = "red", linetype = "dashed", linewidth = 1) +
    labs(title = paste("PBS Distribution -", design_name),
         x = "PBS Value", y = "Frequency") +
    theme_minimal()
  
  # PBS vs position relationship (example: first chromosome)
  chr1_data <- pbs_data %>% filter(chr == unique(pbs_data$chr)[1])
  if (nrow(chr1_data) > 0) {
    p2 <- ggplot(chr1_data, aes(x = start, y = PBS)) +
      geom_point(alpha = 0.5, size = 0.5) +
      geom_smooth(method = "loess", span = 0.1, color = "red") +
      labs(title = paste("PBS Distribution Along Chromosome -", design_name, "-", unique(chr1_data$chr)[1]),
           x = "Position", y = "PBS Value") +
      theme_minimal()
    
    # Save combined plot
    combined_plot <- cowplot::plot_grid(p1, p2, ncol = 1)
    ggsave(file.path(plot_dir, paste0("PBS_plots_", design_name, ".png")), 
           combined_plot, width = 10, height = 8, dpi = 300)
  } else {
    ggsave(file.path(plot_dir, paste0("PBS_distribution_", design_name, ".png")), 
           p1, width = 8, height = 6, dpi = 300)
  }
  
  cat("Visualization plots saved to:", plot_dir, "\n")
}

# Uncomment the following code if visualization is needed
# if (!is.null(analysis_results$pbs_results)) {
#   cat("\nGenerating visualization plots...\n")
#   for (design_name in names(analysis_results$pbs_results)) {
#     plot_PBS_distribution(analysis_results$pbs_results[[design_name]], 
#                           design_name, output_dir)
#   }
# }

cat("\n" %+% "="*60 %+% "\n")
cat("Script execution completed!\n")
cat("Result files saved in:", output_dir, "\n")
cat("="*60, "\n")
