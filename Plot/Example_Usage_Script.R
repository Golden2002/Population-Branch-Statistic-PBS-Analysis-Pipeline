#### Execution Section
# Method 1: Directly run the entire script
source("path/to/main/PBS_visualization_script.R")

# Method 2: Open script file in R Studio, then press Ctrl+Shift+Enter to run entire script

########----Basic Usage----########
# Example 1: Run complete analysis
results <- run_pbs_visualization(
  design_name = "Target_Ref1_Ref2",
  results_dir = "your_results_directory",
  output_dir = "PBS_visualizations"
)

# View results
names(results)
results$stats
results$plots$manhattan
results$plots$distribution

# Example 2: Analyze specific target position
results_with_position <- run_pbs_visualization(
  design_name = "Target_Ref1_Ref2",
  target_position = "chr1:1000000",
  flanking_size = 100000,
  results_dir = "your_results_directory"
)

# Example 3: Quick visualization
quick_results <- quick_visualize(
  design_name = "Target_Ref1_Ref2",
  target_position = "chr1:1000000"
)

########----Advanced Usage----########
# Configuration
design="Sample_Design_1"
dir="path/to/PBS_analysis_results"

# 1. Load data
pbs_data <- load_pbs_results(design, dir)
top_data <- load_top_results(design, dir)

# 2. View statistical information
stats <- calculate_pbs_statistics(pbs_data, threshold_percentile = 0.001)
print(stats)

# 3. Create regional plot
regional_plot <- plot_regional_pbs(
  pbs_data = pbs_data,
  chrom = "chr1",
  center_position = 1715001,
  flanking_size = 100000,
  design_name = design,
  output_dir = "custom_plots"
)

# 4. Create Manhattan plot
manhattan_plot <- plot_pbs_manhattan(
  pbs_data = pbs_data,
  top_data = top_data,
  design_name = design,
  threshold_percentile = 0.001,
  output_dir = "custom_plots"
)

# 5. Create distribution plot
distribution_plot <- plot_pbs_distribution(
  pbs_data = pbs_data,
  design_name = "Target_Ref1_Ref2",
  threshold_percentile = 0.01,
  output_dir = "custom_plots"
)

#######----Diagnostic Methods----########
# 1. Check if files exist
list.files("your_results_directory", pattern = ".*\\.csv$")

# 2. View PBS data structure
pbs_data <- load_pbs_results("Target_Ref1_Ref2", "your_results_directory")
str(pbs_data)
head(pbs_data)

# 3. Check column names
colnames(pbs_data)

# 4. Validate position format
parse_genomic_position("chr1:1000000")
parse_genomic_position("chr2:500000-600000")


########--------########
# Interactive Exploration
########--------########

# 1. Run example to understand functionality
example_results <- example_usage()

# 2. Explore different parameters
# Change threshold
results_5pct <- run_pbs_visualization(
  design_name = "Target_Ref1_Ref2",
  threshold_percentile = 0.05,
  output_dir = "5percent_threshold"
)

# Change region size
results_large_region <- run_pbs_visualization(
  design_name = "Target_Ref1_Ref2",
  target_position = "chr1:1000000",
  flanking_size = 200000,
  output_dir = "large_region"
)

# Change plot parameters
results_custom_plot <- run_pbs_visualization(
  design_name = "Target_Ref1_Ref2",
  point_size = 2,
  text_size = 14,
  output_dir = "custom_plots"
)

########--------########
# Advanced Features
########--------########

# 1. Batch process multiple designs
designs <- c("Design1", "Design2", "Design3")
all_results <- list()

for (design in designs) {
  message("Processing design: ", design)
  results <- run_pbs_visualization(
    design_name = design,
    output_dir = paste0("Visualizations_", design),
    show_plots = FALSE  # Don't display each plot, just save
  )
  all_results[[design]] <- results
}

# 2. Compare PBS distributions across designs
comparison_data <- data.frame()
for (design_name in names(all_results)) {
  design_data <- all_results[[design_name]]$pbs_data
  design_data$design <- design_name
  comparison_data <- rbind(comparison_data, design_data)
}

# Create comparison plot
ggplot(comparison_data, aes(x = PBS, fill = design)) +
  geom_density(alpha = 0.5) +
  labs(title = "PBS Distribution Comparison",
       x = "PBS Value",
       y = "Density") +
  theme_minimal()

# 3. Extract and save top regions
top_regions <- data.frame()
for (design_name in names(all_results)) {
  design_top <- all_results[[design_name]]$top_data
  if (!is.null(design_top)) {
    design_top$design <- design_name
    top_regions <- rbind(top_regions, design_top)
  }
}

# Save top regions
fwrite(top_regions, "all_top_regions.csv")

########--------########
# New Feature: Global Configuration
########--------########


# 1. Set global configuration
config <- setup_configuration(
  results_dir = "my_pbs_results",
  vis_dir = "my_visualizations",
  design_name = "My_Design",
  flanking_size = 100000,
  threshold_percentile = 0.05,
  point_size = 2,
  text_size = 14
)

# 2. Run analysis using configuration
results <- run_pbs_visualization(config = config)

# 3. Modify and re-run
config$flanking_size <- 200000
config$threshold_percentile <- 0.01
results_updated <- run_pbs_visualization(config = config)