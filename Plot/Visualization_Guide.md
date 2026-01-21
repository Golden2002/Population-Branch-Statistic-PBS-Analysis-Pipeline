# PBS Visualization Toolkit for Population Genomics

## Overview

This repository provides a comprehensive R-based toolkit for visualizing and analyzing Population Branch Statistic (PBS) results derived from Fst calculations. The toolkit enables researchers to identify and visualize genomic regions under selection through user-friendly functions that generate publication-quality visualizations.

## Script Architecture

### Main Script: `PBS_Visualization.R`
Contains the core functionality including:
- Data loading and preprocessing functions
- Statistical calculation routines
- Multiple visualization options (Manhattan plots, regional views, distribution plots)
- Configuration management system
- Batch processing capabilities

### Example Script: `Example_Usage.R`
Provides practical examples demonstrating:
- Basic and advanced usage patterns
- Parameter exploration workflows
- Batch processing templates
- Interactive exploration methods

## Requirements

### R Packages
```r
ggplot2      # Data visualization
dplyr        # Data manipulation
data.table   # Fast data loading
ggrepel      # Intelligent text labeling
patchwork    # Multi-panel plot arrangement
viridis      # Colorblind-friendly palettes
gridExtra    # Advanced plot arrangement
```

Install required packages:
```r
install.packages(c("ggplot2", "dplyr", "data.table", "ggrepel", 
                   "patchwork", "viridis", "gridExtra"))
```

## Input Data Format

### Required PBS Results File
The toolkit expects PBS calculation results in CSV format with the following structure:

**Mandatory Columns:**
- `chr`: Chromosome identifier (e.g., "chr1", "chr2", "chrX")
- `start`: Start position of genomic window
- `end`: End position of genomic window
- `PBS`: Calculated PBS value

**Example File Structure:**
```
chr,start,end,PBS,fst_TR1,fst_TR2,fst_R1R2
chr1,1000,2000,0.152,0.120,0.135,0.087
chr1,2001,3000,0.087,0.065,0.072,0.045
```

### File Naming Conventions
The toolkit automatically searches for files with these patterns:
- `PBS_full_[design_name].csv`
- `PBS_results_[design_name].csv`
- `[design_name]_PBS.csv`

## Quick Start Guide

### 1. Basic Analysis
```r
# Source the main script
source("PBS_Visualization.R")

# Run complete analysis
results <- run_pbs_visualization(
  design_name = "Target_Ref1_Ref2",
  results_dir = "path/to/pbs/results",
  output_dir = "Visualization_Output"
)
```

### 2. Regional Analysis
```r
# Analyze specific genomic region
results <- run_pbs_visualization(
  design_name = "Target_Ref1_Ref2",
  target_position = "chr1:1715001",
  flanking_size = 100000,  # 100kb flanking region
  threshold_percentile = 0.001  # Top 0.1%
)
```

### 3. Quick Visualization
```r
# Generate visualizations with default settings
quick_results <- quick_visualize(
  design_name = "Jino_Tibetan_Blang",
  target_position = "chr1:1000000"
)
```

## Core Functions

### Data Loading Functions
- `load_pbs_results()`: Load PBS calculation results
- `load_top_results()`: Load top percentile regions
- `parse_genomic_position()`: Parse genomic coordinates

### Visualization Functions
- `plot_regional_pbs()`: Generate regional PBS plots
- `plot_pbs_manhattan()`: Create genome-wide Manhattan plots
- `plot_pbs_distribution()`: Plot PBS value distributions

### Analysis Functions
- `calculate_pbs_statistics()`: Compute summary statistics
- `create_summary_visualization()`: Generate comprehensive visualization set

### Workflow Functions
- `run_pbs_visualization()`: Complete analysis workflow
- `quick_visualize()`: Simplified interface
- `example_usage()`: Run demonstration analysis

## Output Specifications

### Generated Files
1. **Manhattan Plots** (`PBS_manhattan_[design].png`):
   - Genome-wide PBS distribution
   - Top percentile threshold line
   - Chromosome-specific coloring
   - Zoomed versions for high-PBS regions

2. **Regional Plots** (`PBS_region_[position]_[design].png`):
   - PBS signal around target position
   - Target position indicator
   - 99th percentile reference line
   - Position-specific annotations

3. **Distribution Plots** (`PBS_distribution_[design].png`):
   - Histogram of PBS values
   - 99th percentile and threshold lines
   - Summary statistics annotations

4. **Summary Reports**:
   - Statistical summary in console output
   - Configuration settings
   - Analysis metadata

### Return Object Structure
```r
results <- run_pbs_visualization(...)

# Access components
results$config      # Analysis configuration
results$pbs_data    # Loaded PBS data
results$top_data    # Top percentile regions
results$stats       # Statistical summary
results$plots       # Generated plot objects
```

## Advanced Features

### Batch Processing
```r
# Analyze multiple designs
designs <- c("Design1", "Design2", "Design3")
all_results <- list()

for (design in designs) {
  results <- run_pbs_visualization(
    design_name = design,
    output_dir = paste0("Visualizations_", design),
    show_plots = FALSE
  )
  all_results[[design]] <- results
}
```

### Configuration Management
```r
# Set global configuration
config <- setup_configuration(
  results_dir = "my_pbs_results",
  vis_dir = "my_visualizations",
  design_name = "My_Design",
  flanking_size = 200000,
  threshold_percentile = 0.05,
  point_size = 2,
  text_size = 14
)

# Run with configuration
results <- run_pbs_visualization(config = config)
```

### Parameter Exploration
```r
# Compare different thresholds
results_1pct <- run_pbs_visualization(
  design_name = "Target_Ref1_Ref2",
  threshold_percentile = 0.01,
  output_dir = "1percent_threshold"
)

results_5pct <- run_pbs_visualization(
  design_name = "Target_Ref1_Ref2",
  threshold_percentile = 0.05,
  output_dir = "5percent_threshold"
)
```

## Troubleshooting Guide

### Common Issues

1. **File Not Found Errors**
   ```r
   # Check available files
   list.files("your_results_directory", pattern = ".*\\.csv$")
   
   # Verify design name matches file naming
   available_files <- list.files(results_dir, pattern = ".*\\.csv$")
   ```

2. **Data Structure Issues**
   ```r
   # Examine data structure
   pbs_data <- load_pbs_results(design_name, results_dir)
   str(pbs_data)
   head(pbs_data)
   
   # Verify column names
   colnames(pbs_data)
   ```

3. **Visualization Problems**
   - Ensure `midpoint` column exists: `pbs_data$midpoint <- (pbs_data$start + pbs_data$end)/2`
   - Check chromosome naming format: Should be "chr1", "chr2", etc.
   - Verify PBS values are numeric: `class(pbs_data$PBS)`

### Position Format Validation
```r
# Valid formats
parse_genomic_position("chr1:1000000")
parse_genomic_position("chr2:500000-600000")

# Invalid formats
parse_genomic_position("rs12345")  # RS IDs not supported
parse_genomic_position("1:1000000")  # Missing 'chr' prefix
```

## Best Practices

### 1. Data Preparation
- Ensure PBS results include required columns
- Use consistent chromosome naming ("chr1", not "1" or "Chr1")
- Remove missing values before analysis

### 2. Parameter Selection
- Use `threshold_percentile = 0.001` for stringent selection
- Adjust `flanking_size` based on genomic context
- Consider `point_size` and `text_size` for publication figures

### 3. Output Management
- Specify unique `output_dir` for different analyses
- Use descriptive `design_name` for easy identification
- Save summary statistics for reproducibility

### 4. Quality Control
```r
# Verify data quality
stats <- calculate_pbs_statistics(pbs_data)
print(stats)

# Check for outliers
boxplot(pbs_data$PBS, main = "PBS Value Distribution")
```

## Citation and Attribution

When using this toolkit in publications, please cite:

> Population Branch Statistic (PBS) Visualization Toolkit. Available at: [GitHub Repository URL]

## Support and Contributions

For issues, questions, or contributions:
1. Report bugs via GitHub Issues
2. Submit pull requests for enhancements
3. Contact maintainers for technical support

## License

This project is licensed under the MIT License - see the LICENSE file for details.

---

**Note**: This toolkit is designed for research use. Always validate results with appropriate statistical methods and consider biological context in interpretation.
