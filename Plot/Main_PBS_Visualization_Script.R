# ============================================================================
# PBS Visualization Script for R Studio (Corrected Version)
# ============================================================================

library(ggplot2)
library(dplyr)
library(data.table)
library(ggrepel)
library(patchwork)
library(viridis)
library(gridExtra)

# ============================================================================
# 1. Global Configuration Functions
# ============================================================================

setup_configuration <- function(results_dir = ".",
                                vis_dir = "PBS_visualizations",
                                design_name = "Target_Ref1_Ref2",
                                flanking_size = 50000,
                                threshold_percentile = 0.01,
                                point_size = 1.5,
                                text_size = 12) {
  
  if (!dir.exists(vis_dir)) {
    dir.create(vis_dir, recursive = TRUE)
    message("Created output directory: ", vis_dir)
  }
  
  config <- list(
    results_dir = results_dir,
    vis_dir = vis_dir,
    design_name = design_name,
    flanking_size = flanking_size,
    threshold_percentile = threshold_percentile,
    point_size = point_size,
    text_size = text_size
  )
  
  return(config)
}

# ============================================================================
# 2. Data Loading Functions
# ============================================================================

load_pbs_results <- function(design_name, results_dir = ".") {
  
  # Find files
  available_files <- list.files(results_dir, pattern = ".*\\.csv$", full.names = FALSE)
  message("CSV files found in ", results_dir, ": ", 
          paste(available_files, collapse = ", "))
  
  # Possible filename patterns
  possible_files <- c(
    paste0("PBS_full_", design_name, ".csv"),
    paste0("PBS_results_", design_name, ".csv"),
    paste0(design_name, "_PBS.csv")
  )
  
  pbs_file <- NULL
  for (file in possible_files) {
    file_path <- file.path(results_dir, file)
    if (file.exists(file_path)) {
      pbs_file <- file_path
      break
    }
  }
  
  if (is.null(pbs_file)) {
    # Try to find any file containing design name
    all_files <- list.files(results_dir, pattern = ".*\\.csv$", full.names = TRUE)
    matching_files <- all_files[grepl(design_name, all_files)]
    if (length(matching_files) > 0) {
      pbs_file <- matching_files[1]
    } else {
      stop(paste("PBS results file for design", design_name, "not found"))
    }
  }
  
  message("Loading PBS results from: ", pbs_file)
  pbs_data <- fread(pbs_file)
  
  # Basic data information
  message("Data loaded successfully.")
  message("Number of rows: ", nrow(pbs_data))
  message("Number of chromosomes: ", length(unique(pbs_data$chr)))
  message("Column names: ", paste(names(pbs_data), collapse = ", "))
  
  # Calculate window center point (for plotting)
  if ("start" %in% names(pbs_data) && "end" %in% names(pbs_data)) {
    pbs_data <- pbs_data %>%
      mutate(midpoint = (start + end) / 2)
    message("Added midpoint column for plotting.")
  }
  
  return(pbs_data)
}

load_top_results <- function(design_name, results_dir = ".") {
  
  # Find top results file
  possible_files <- list.files(results_dir, 
                               pattern = paste0(".*top.*", design_name, ".*\\.csv$"), 
                               full.names = TRUE)
  
  if (length(possible_files) > 0) {
    top_file <- possible_files[1]
    message("Loading top results from: ", top_file)
    top_data <- fread(top_file)
    
    # Calculate window center point (if relevant columns exist)
    if ("start" %in% names(top_data) && "end" %in% names(top_data)) {
      top_data <- top_data %>%
        mutate(midpoint = (start + end) / 2)
    }
    
    message("Loaded ", nrow(top_data), " top percentile windows.")
    return(top_data)
  } else {
    message("Warning: Top results file for design ", design_name, " not found")
    return(NULL)
  }
}

parse_genomic_position <- function(position_string) {
  
  # Remove whitespace
  position_string <- gsub("\\s", "", position_string)
  
  # Check format
  if (grepl(":", position_string)) {
    parts <- strsplit(position_string, ":")[[1]]
    chrom <- parts[1]
    
    if (grepl("-", parts[2])) {
      range_parts <- strsplit(parts[2], "-")[[1]]
      start_pos <- as.numeric(range_parts[1])
      end_pos <- as.numeric(range_parts[2])
    } else {
      start_pos <- as.numeric(parts[2])
      end_pos <- start_pos  # Single position
    }
    
    message("Parsed position: ", chrom, ":", start_pos, "-", end_pos)
    return(list(chrom = chrom, start = start_pos, end = end_pos))
    
  } else if (grepl("^rs", position_string)) {
    stop("RS ID conversion not implemented. Please provide genomic position in format 'chr:position'.")
  } else {
    stop("Invalid position format. Use 'chr:position' or 'chr:start-end'")
  }
}

calculate_pbs_statistics <- function(pbs_data, threshold_percentile = 0.01) {
  
  stats <- list(
    n_windows = nrow(pbs_data),
    n_chromosomes = length(unique(pbs_data$chr)),
    mean_pbs = mean(pbs_data$PBS, na.rm = TRUE),
    median_pbs = median(pbs_data$PBS, na.rm = TRUE),
    min_pbs = min(pbs_data$PBS, na.rm = TRUE),
    max_pbs = max(pbs_data$PBS, na.rm = TRUE),
    sd_pbs = sd(pbs_data$PBS, na.rm = TRUE),
    threshold_value = quantile(pbs_data$PBS, 1 - threshold_percentile, na.rm = TRUE),
    n_above_threshold = sum(pbs_data$PBS > quantile(pbs_data$PBS, 1 - threshold_percentile, na.rm = TRUE), na.rm = TRUE)
  )
  
  return(stats)
}

# ============================================================================
# 3. Visualization Functions (Corrected)
# ============================================================================

#' Prepare Manhattan plot data
#' 
#' @param pbs_data PBS data frame
#' @return Processed data frame

prepare_manhattan_data <- function(pbs_data) {
  
  # Sort by chromosome
  # Supports chr1, chr2 ... chrX, chrY, chrM, chrMT
  chrom_order <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM", "chrMT")
  
  # Convert chromosome to factor
  pbs_data$chr <- factor(pbs_data$chr, levels = chrom_order)
  
  # Remove chromosomes not in standard order
  pbs_data <- pbs_data %>% filter(!is.na(chr))
  
  # Calculate chromosome index
  pbs_data <- pbs_data %>%
    mutate(chr_index = as.numeric(chr))  # Directly use factor numeric values
  
  if (all(is.na(pbs_data$chr_index))) {
    stop("chr_index is all NA. Please check chromosome naming follows 'chr1','chr2',...'chrX' format")
  }
  
  
  # Calculate chromosome center positions (for x-axis labels)
  axis_df <- pbs_data %>%
    group_by(chr) %>%
    summarize(
      center = mean(chr_index, na.rm = TRUE)  # Use mean of chr_index
    ) %>%
    ungroup()
  
  return(list(data = pbs_data, axis_df = axis_df))
}

plot_regional_pbs <- function(pbs_data, chrom, center_position, flanking_size = 50000, 
                              design_name = "Unknown", output_dir = ".",
                              point_size = 1.5, text_size = 12, show_plot = TRUE) {
  
  # Define region boundaries
  region_start <- max(0, center_position - flanking_size)
  region_end <- center_position + flanking_size
  
  message("Plotting region: ", chrom, ":", region_start, "-", region_end)
  
  # Filter data for specified region
  # Note: Ensure chrom column matches input chromosome name
  if (is.factor(pbs_data$chr)) {
    region_data <- pbs_data %>%
      filter(as.character(chr) == chrom & 
               midpoint >= region_start & 
               midpoint <= region_end)
  } else {
    region_data <- pbs_data %>%
      filter(chr == chrom & 
               midpoint >= region_start & 
               midpoint <= region_end)
  }
  
  if (nrow(region_data) == 0) {
    warning("No PBS data found in specified region.")
    return(NULL)
  }
  
  # Calculate PBS statistics for this region
  region_stats <- calculate_pbs_statistics(region_data)
  
  # Create plot
  p <- ggplot(region_data, aes(x = midpoint, y = PBS)) +
    geom_line(color = "steelblue", linewidth = 1) +
    geom_point(color = "steelblue", size = point_size, alpha = 0.7) +
    geom_vline(xintercept = center_position, color = "red", 
               linetype = "dashed", linewidth = 1, alpha = 0.7) +
    geom_hline(yintercept = quantile(pbs_data$PBS, 0.99, na.rm = TRUE), 
               color = "darkorange", linetype = "dotted", linewidth = 0.8) +
    labs(
      title = paste("PBS Signal - ", chrom, ":", format(center_position, big.mark = ",")),
      subtitle = paste("Design: ", design_name, 
                       " | Region: ", format(region_start, big.mark = ","), "-", 
                       format(region_end, big.mark = ",")),
      x = "Genomic Position",
      y = "PBS Value",
      caption = paste("Max PBS in region: ", round(region_stats$max_pbs, 4),
                      " | Mean: ", round(region_stats$mean_pbs, 4))
    ) +
    theme_minimal(base_size = text_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "gray50"),
      plot.caption = element_text(hjust = 0.5, color = "gray40", size = text_size * 0.8),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
    ) +
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    annotate("text", 
             x = center_position, 
             y = max(region_data$PBS, na.rm = TRUE) * 0.95,
             label = "Target Position",
             color = "red",
             hjust = -0.1,
             vjust = 0)
  
  # Save image
  if (!is.null(output_dir)) {
    filename <- paste0("PBS_region_", gsub("[:]", "_", chrom), "_", 
                       center_position, "_", design_name, ".png")
    filepath <- file.path(output_dir, filename)
    ggsave(filepath, p, width = 12, height = 6, dpi = 300, bg = "white")
    message("Regional plot saved to: ", filepath)
  }
  
  # Display image (if requested)
  if (show_plot) {
    print(p)
  }
  
  return(p)
}

plot_pbs_manhattan <- function(pbs_data, top_data = NULL, design_name = "Unknown", 
                               threshold_percentile = 0.01, output_dir = ".",
                               point_size = 1.5, text_size = 12, show_plot = TRUE) {
  
  # Prepare Manhattan plot data
  message("Preparing Manhattan plot data...")
  manhattan_data <- prepare_manhattan_data(pbs_data)
  pbs_data_processed <- manhattan_data$data
  axis_df <- manhattan_data$axis_df
  
  # Calculate PBS threshold (for highlighting)
  pbs_threshold <- quantile(pbs_data$PBS, 1 - threshold_percentile, na.rm = TRUE)
  
  message("Manhattan plot threshold (top ", threshold_percentile*100, "%): ", 
          round(pbs_threshold, 4))
  
  # Create Manhattan plot
  p <- ggplot(pbs_data_processed, aes(x = chr_index, y = PBS)) +
    # Add alternating chromosome background
    geom_rect(data = axis_df,
              aes(xmin = as.numeric(chr) - 0.5, 
                  xmax = as.numeric(chr) + 0.5,
                  ymin = -Inf, ymax = Inf,
                  fill = ifelse(as.numeric(chr) %% 2 == 0, "even", "odd")),
              alpha = 0.2, inherit.aes = FALSE) + # Add inherit.aes = FALSE
    scale_fill_manual(values = c("even" = "gray90", "odd" = "white"), guide = "none") +
    
    # Add all points
    geom_point(aes(color = chr), alpha = 0.6, size = point_size, position = position_jitter(width = 0.3)) +
    
    # Add threshold line
    geom_hline(yintercept = pbs_threshold, 
               color = "red", 
               linetype = "dashed", 
               linewidth = 0.8,
               alpha = 0.7) +
    
    # Customize plot appearance
    labs(
      title = paste("PBS Manhattan Plot - ", design_name),
      subtitle = paste("Highlighting top ", threshold_percentile*100, "% PBS values (threshold = ",
                       round(pbs_threshold, 4), ")", sep = ""),
      x = "Chromosome",
      y = "PBS Value"
    ) +
    theme_minimal(base_size = text_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "gray50"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    
    # Set x-axis and chromosome labels
    scale_x_continuous(
      breaks = axis_df$center,
      labels = axis_df$chr,
      expand = expansion(mult = c(0.01, 0.01)),
      limits = c(0.5, max(axis_df$center) + 0.5)  # Add explicit limits
    ) +
    
    # Set point colors by chromosome
    scale_color_viridis(discrete = TRUE, option = "plasma") +
    
    # Add threshold annotation
    annotate("text",
             x = max(axis_df$center, na.rm = TRUE) * 0.98,  # Use axis_df$center instead of pbs_data_processed$chr_index
             y = pbs_threshold * 1.02,
             label = paste("Top", threshold_percentile*100, "% threshold"),
             color = "red",
             size = 3.5,
             hjust = 1)
  
  # If top data provided, highlight
  if (!is.null(top_data) && nrow(top_data) > 0) {
    message("Highlighting top regions...")
    
    # Prepare top data
    top_data_processed <- prepare_manhattan_data(top_data)$data
    
    # Add highlight points
    p <- p + 
      geom_point(data = top_data_processed, 
                 aes(x = chr_index, y = PBS),
                 color = "red", 
                 size = point_size * 1.5,
                 alpha = 0.8)
    
    # Label some top points
    if (nrow(top_data_processed) > 0) {
      # Select top 10 regions by PBS for labeling
      top_to_label <- top_data_processed %>%
        arrange(desc(PBS)) %>%
        head(min(10, nrow(top_data_processed))) %>%
        mutate(label = paste(as.character(chr), ":", 
                             format(round(midpoint), big.mark = ","), sep = ""))
      
      if (nrow(top_to_label) > 0) {
        p <- p + 
          geom_label_repel(data = top_to_label,
                           aes(x = chr_index, y = PBS, label = label),
                           size = 3,
                           box.padding = 0.5,
                           point.padding = 0.3,
                           segment.color = "gray50",
                           min.segment.length = 0,
                           max.overlaps = Inf)
      }
    }
  }
  
  # Save image
  if (!is.null(output_dir)) {
    filename <- paste0("PBS_manhattan_", design_name, ".png")
    filepath <- file.path(output_dir, filename)
    ggsave(filepath, p, width = 16, height = 8, dpi = 300, bg = "white")
    message("Manhattan plot saved to: ", filepath)
  }
  
  # Display image (if requested)
  if (show_plot) {
    print(p)
  }
  
  # Create zoomed version (focused on top signals)
  if (!is.null(top_data) && nrow(top_data) > 0) {
    max_pbs <- max(top_data$PBS, na.rm = TRUE)
    p_zoom <- p + 
      coord_cartesian(ylim = c(pbs_threshold * 0.9, max_pbs * 1.1)) +
      labs(title = paste("PBS Manhattan Plot (Zoomed) - ", design_name),
           subtitle = paste("Focusing on top", threshold_percentile*100, "% signals"))
    
    if (!is.null(output_dir)) {
      filename_zoom <- paste0("PBS_manhattan_zoomed_", design_name, ".png")
      filepath_zoom <- file.path(output_dir, filename_zoom)
      ggsave(filepath_zoom, p_zoom, width = 16, height = 8, dpi = 300, bg = "white")
      message("Zoomed Manhattan plot saved to: ", filepath_zoom)
    }
    
    if (show_plot) {
      print(p_zoom)
    }
  }
  
  return(p)
}

plot_pbs_distribution <- function(pbs_data, design_name = "Unknown", 
                                  threshold_percentile = 0.01, output_dir = ".",
                                  text_size = 12, show_plot = TRUE) {
  
  # Calculate thresholds
  pbs_threshold <- quantile(pbs_data$PBS, 1 - threshold_percentile, na.rm = TRUE)
  pbs_99th <- quantile(pbs_data$PBS, 0.99, na.rm = TRUE)
  
  # Create distribution plot
  dist_plot <- ggplot(pbs_data, aes(x = PBS)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = pbs_99th, 
               color = "red", linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = pbs_threshold, 
               color = "darkorange", linetype = "dotted", linewidth = 1) +
    labs(
      title = paste("PBS Value Distribution - ", design_name),
      subtitle = paste("Mean = ", round(mean(pbs_data$PBS, na.rm = TRUE), 4),
                       ", Median = ", round(median(pbs_data$PBS, na.rm = TRUE), 4)),
      x = "PBS Value",
      y = "Frequency"
    ) +
    theme_minimal(base_size = text_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "gray50")
    ) +
    annotate("text",
             x = pbs_99th,
             y = Inf,
             label = "99th percentile",
             color = "red",
             hjust = -0.1,
             vjust = 2) +
    annotate("text",
             x = pbs_threshold,
             y = Inf,
             label = paste("Top", threshold_percentile*100, "% (", 
                           round(pbs_threshold, 4), ")"),
             color = "darkorange",
             hjust = -0.1,
             vjust = 4)
  
  # Save image
  if (!is.null(output_dir)) {
    filename_dist <- paste0("PBS_distribution_", design_name, ".png")
    filepath_dist <- file.path(output_dir, filename_dist)
    ggsave(filepath_dist, dist_plot, width = 10, height = 6, dpi = 300, bg = "white")
    message("Distribution plot saved to: ", filepath_dist)
  }
  
  # Display image (if requested)
  if (show_plot) {
    print(dist_plot)
  }
  
  return(dist_plot)
}

create_summary_visualization <- function(pbs_data, top_data = NULL, design_name = "Unknown", 
                                         target_position = NULL, config, show_plots = TRUE) {
  
  plots <- list()
  
  # Extract parameters from configuration
  output_dir <- config$vis_dir
  flanking_size <- config$flanking_size
  threshold_percentile <- config$threshold_percentile
  point_size <- config$point_size
  text_size <- config$text_size
  
  # Create Manhattan plot
  message("Creating Manhattan plot...")
  manhattan_plot <- plot_pbs_manhattan(pbs_data, top_data, design_name, 
                                       threshold_percentile, output_dir,
                                       point_size, text_size, show_plots)
  plots$manhattan <- manhattan_plot
  
  # If target position provided, create regional plot
  if (!is.null(target_position)) {
    message("Creating regional plot...")
    pos_info <- parse_genomic_position(target_position)
    regional_plot <- plot_regional_pbs(pbs_data, 
                                       pos_info$chrom, 
                                       pos_info$start,
                                       flanking_size,
                                       design_name,
                                       output_dir,
                                       point_size,
                                       text_size,
                                       show_plots)
    plots$regional <- regional_plot
    
    # Create combined plot
    if (!is.null(regional_plot)) {
      combined_plot <- regional_plot / manhattan_plot + 
        plot_layout(heights = c(1, 2)) +
        plot_annotation(
          title = paste("PBS Analysis Summary - ", design_name),
          subtitle = paste("Regional view for target position ", target_position, " and genome-wide Manhattan plot"),
          tag_levels = 'A'
        )
      
      if (!is.null(output_dir)) {
        filename <- paste0("PBS_summary_", design_name, "_", 
                           gsub("[^a-zA-Z0-9]", "_", target_position), ".png")
        filepath <- file.path(output_dir, filename)
        ggsave(filepath, combined_plot, width = 14, height = 10, dpi = 300, bg = "white")
        message("Summary visualization saved to: ", filepath)
      }
      
      if (show_plots) {
        print(combined_plot)
      }
      
      plots$combined <- combined_plot
    }
  }
  
  # Create PBS distribution plot
  message("Creating PBS distribution plot...")
  dist_plot <- plot_pbs_distribution(pbs_data, design_name, 
                                     threshold_percentile, output_dir,
                                     text_size, show_plots)
  plots$distribution <- dist_plot
  
  return(plots)
}

# ============================================================================
# 4. Interactive Analysis Functions
# ============================================================================

run_pbs_visualization <- function(results_dir = ".",
                                  design_name = "Target_Ref1_Ref2",
                                  target_position = NULL,
                                  flanking_size = 50000,
                                  threshold_percentile = 0.01,
                                  output_dir = "PBS_visualizations",
                                  point_size = 1.5,
                                  text_size = 12,
                                  show_plots = TRUE) {
  
  message(strrep("=", 60))
  message("PBS Visualization Script (R Studio Version)")
  message(strrep("=", 60))
  message("")
  
  # Setup configuration
  message("Setting up configuration...")
  config <- setup_configuration(
    results_dir = results_dir,
    vis_dir = output_dir,
    design_name = design_name,
    flanking_size = flanking_size,
    threshold_percentile = threshold_percentile,
    point_size = point_size,
    text_size = text_size
  )
  
  # Display configuration
  message("Configuration:")
  message("  Results directory: ", config$results_dir)
  message("  Output directory: ", config$vis_dir)
  message("  Design: ", config$design_name)
  message("  Region size: ", config$flanking_size, " bp")
  message("  Threshold percentile: ", config$threshold_percentile)
  if (!is.null(target_position)) {
    message("  Target position: ", target_position)
  }
  message("")
  
  # Load data
  message("Loading PBS results...")
  pbs_data <- load_pbs_results(config$design_name, config$results_dir)
  
  # Load top results
  message("Loading top percentile results...")
  top_data <- load_top_results(config$design_name, config$results_dir)
  
  # Calculate and display basic statistics
  stats <- calculate_pbs_statistics(pbs_data, config$threshold_percentile)
  message("")
  message("PBS Statistics:")
  message("  Number of windows: ", stats$n_windows)
  message("  Number of chromosomes: ", stats$n_chromosomes)
  message("  Mean PBS: ", round(stats$mean_pbs, 6))
  message("  Median PBS: ", round(stats$median_pbs, 6))
  message("  Minimum PBS: ", round(stats$min_pbs, 6))
  message("  Maximum PBS: ", round(stats$max_pbs, 6))
  message("  Standard deviation: ", round(stats$sd_pbs, 6))
  message("  Top ", config$threshold_percentile*100, "% threshold: ", round(stats$threshold_value, 6))
  message("  Windows above threshold: ", stats$n_above_threshold)
  message("")
  
  # Create visualizations
  message("Creating visualizations...")
  plots <- create_summary_visualization(
    pbs_data = pbs_data,
    top_data = top_data,
    design_name = config$design_name,
    target_position = target_position,
    config = config,
    show_plots = show_plots
  )
  
  message("")
  message(strrep("=", 60))
  message("Visualization completed!")
  message("All images saved to: ", config$vis_dir)
  message(strrep("=", 60))
  
  # Return results for further analysis
  results <- list(
    config = config,
    pbs_data = pbs_data,
    top_data = top_data,
    stats = stats,
    plots = plots
  )
  
  return(results)
}

quick_visualize <- function(design_name, target_position = NULL, results_dir = ".") {
  message("Quick visualization for design: ", design_name)
  
  results <- run_pbs_visualization(
    results_dir = results_dir,
    design_name = design_name,
    target_position = target_position,
    output_dir = paste0("Quick_", design_name, "_Visualizations"),
    show_plots = TRUE
  )
  
  return(results)
}

# ============================================================================
# 5. Example Usage Functions
# ============================================================================

example_usage <- function() {
  
  message("Example 1: Basic visualization for a design")
  message("--------------------------------------------")
  message("results1 <- run_pbs_visualization(")
  message('  design_name = "Target_Ref1_Ref2",')
  message('  results_dir = ".",')
  message('  output_dir = "PBS_visualizations"')
  message(")")
  message("")
  
  message("Example 2: Visualization with target position")
  message("---------------------------------------------")
  message("results2 <- run_pbs_visualization(")
  message('  design_name = "Target_Ref1_Ref2",')
  message('  target_position = "chr1:1000000",')
  message('  flanking_size = 100000,')
  message('  threshold_percentile = 0.05')
  message(")")
  message("")
  
  message("Example 3: Quick visualization")
  message("-------------------------------")
  message('results3 <- quick_visualize("Target_Ref1_Ref2", "chr1:1000000")')
  message("")
  
  message("Example 4: Custom visualization")
  message("------------------------------------")
  message("results4 <- run_pbs_visualization(")
  message('  design_name = "Target_Ref1_Ref2",')
  message('  target_position = "chr1:1000000",')
  message('  flanking_size = 50000,')
  message('  threshold_percentile = 0.01,')
  message('  point_size = 2,')
  message('  text_size = 14,')
  message('  show_plots = TRUE')
  message(")")
  
  # Actually run an example
  message("")
  message("Running Example 1 (basic visualization)...")
  message("")
  
  # Check if example data exists
  example_file <- paste0("PBS_full_Target_Ref1_Ref2.csv")
  if (!file.exists(example_file)) {
    message("Example file not found. Creating dummy data for demonstration...")
    
    # Create dummy data for demonstration
    set.seed(123)
    dummy_data <- data.frame(
      chr = rep(c("1", "2", "3"), each = 100),
      start = rep(seq(1, 1000000, length.out = 100), 3),
      end = rep(seq(10000, 1010000, length.out = 100), 3),
      PBS = c(rnorm(100, mean = 0.1, sd = 0.05),
              rnorm(100, mean = 0.15, sd = 0.08),
              rnorm(100, mean = 0.08, sd = 0.03))
    )
    
    dummy_data <- dummy_data %>%
      mutate(
        fst_TR1 = PBS * 0.8,
        fst_TR2 = PBS * 0.9,
        fst_R1R2 = PBS * 0.2,
        T_TR1 = -log(1 - fst_TR1),
        T_TR2 = -log(1 - fst_TR2),
        T_R1R2 = -log(1 - fst_R1R2),
        design = "Target_Ref1_Ref2",
        midpoint = (start + end) / 2
      )
    
    # Save dummy data
    fwrite(dummy_data, example_file)
    message("Dummy data saved to: ", example_file)
  }
  
  # Run example
  results <- run_pbs_visualization(
    design_name = "Target_Ref1_Ref2",
    results_dir = ".",
    output_dir = "Example_Visualizations",
    show_plots = TRUE
  )
  
  message("")
  message("Example completed!")
  message("Results stored in 'results' object.")
  message("You can access:")
  message("  - PBS data: results$pbs_data")
  message("  - Statistics: results$stats")
  message("  - Plots: results$plots")
  message("  - Configuration: results$config")
  
  return(results)
}

# ============================================================================
# 6. Main Execution Section
# ============================================================================

# When loading this script in R Studio, you can run:
# 1. example_usage() - See demonstration
# 2. run_pbs_visualization() - Use your parameters
# 3. Use individual functions directly

message("PBS visualization script loaded successfully.")
message("")
message("Available functions:")
message("  1. run_pbs_visualization() - Main function for interactive analysis")
message("  2. quick_visualize() - Simplified interface")
message("  3. example_usage() - Run example analysis")
message("  4. load_pbs_results() - Load PBS data")
message("  5. plot_regional_pbs() - Create regional PBS plot")
message("  6. plot_pbs_manhattan() - Create Manhattan plot")
message("  7. plot_pbs_distribution() - Create distribution plot")
message("")
message("To get started, run: example_usage()")