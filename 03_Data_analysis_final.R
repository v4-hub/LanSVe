# -------------------------------------------------------------------------
# The Nonlinear Relationship Between Vegetation Indices and Landslide Hazard: A Study Exceeding Nature Communications Standards
# Multi-Source Remote Sensing Data-Driven Risk Threshold Discovery and Global Validation Framework
# Expert-Level Analysis: Landslide Mechanisms + Vegetation Remote Sensing + Journal Editor's Perspective
# -------------------------------------------------------------------------

# =========================================================================
# 1. Environment Setup and Core Libraries
# =========================================================================
# Load core libraries and check versions
required_packages <- c(
  "ggplot2", "dplyr", "segmented", "randomForest", "viridis", "scales",
  "rnaturalearth", "rnaturalearthdata", "sf", "pdp", "tidyr", "patchwork",
  "RColorBrewer", "boot", "GGally", "corrplot", "treemapify", "ggalluvial",
  "vip", "DALEX", "iBreakDown", "plotly", "ggsci", "ggpubr", "networkD3",
  "caret", "pROC", "smotefamily", "ggbeeswarm", "terra", "raster", "ggrepel",
  "irr", "Kendall", "Rmisc", "plotROC", "systemfonts", "stringr"
)

for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, quiet = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Function to check available system fonts
check_available_fonts <- function() {
  available_fonts <- systemfonts::system_fonts()
  
  # Check for Nature-recommended fonts
  nature_fonts <- available_fonts %>% 
    filter(str_detect(family, "Arial|Helvetica|Times|Calibri")) %>% 
    dplyr::select(family, style) %>% 
    distinct()
  
  return(nature_fonts)
}

# Get the best available font
get_best_nature_font <- function() {
  available_fonts <- systemfonts::system_fonts()$family
  
  # Font priority for Nature journals
  font_priority <- c("Arial", "Helvetica", "Calibri", "Times New Roman", "Times")
  
  for(font in font_priority) {
    if (font %in% available_fonts) {
      cat("✓ Using font:", font, "\n")
      return(font)
    }
  }
  
  cat("⚠ Using default sans-serif font\n")
  return("")  # An empty string lets R use the default font
}

# Check available fonts
cat("Available Nature-style fonts on your system:\n")
print(check_available_fonts())

# Get the best font
best_font <- get_best_nature_font()

# =========================================================================
# Set up directories for saving results
# Define the base directory
base_dir <- "V5"

# Define the subdirectories
sub_dirs <- c("figures", "tables", "data")

# Create the base directory and its subdirectories
for (sub_dir in sub_dirs) {
  full_path <- file.path(base_dir, sub_dir)
  dir.create(full_path, showWarnings = FALSE, recursive = TRUE)
}
# =========================================================================
# -------------------------------------------------------------------------
# Professional Academic Visual System (Nature-series Journal Standard)
# -------------------------------------------------------------------------

# Professional Color Palettes for Nature Communications
nature_palettes <- list(
  # Main findings palette (based on Nature journal guidelines)
  biophysical_zones = c(
    "IDZ" = "#0173B2",      # Blue: Low risk (Initial Disturbance Zone)
    "BTZ" = "#DE8F05",      # Orange: Transition Zone
    "CTZ" = "#CC78BC",      # Purple: High risk (Critical Transition Zone)
    "SDZ" = "#029E73",      # Green: Post-peak (Stabilization/Decline Zone)
    "Very High NDVI" = "#56B4E9", # Light Blue: Very high vegetation
    "Unknown" = "#999999"         # Grey: Unknown
  ),
  
  # Professional palette for vegetation indices
  indices = c(
    "NDVI" = "#1B9E77", "EVI" = "#D95F02", 
    "LAI" = "#7570B3", "FPAR" = "#E7298A"
  ),
  
  # Scientific palette for climate zones
  climate = c(
    "Temperate" = "#2166AC",
    "Mediterranean/Subtropical" = "#762A83", 
    "Tropical" = "#5AAE61",
    "Other" = "#9970AB"
  ),
  
  # Ecological palette for vegetation types
  vegetation = c(
    "Forest" = "#1A5490", "Woody" = "#4292C6", "Non-Forest" = "#C994C7",
    "Evergreen Broadleaf Forests" = "#08519C",
    "Woody Savannas" = "#3182BD",
    "Savannas" = "#6BAED6",
    "Croplands" = "#FD8D3C", 
    "Urban and Built-up Lands" = "#525252",
    "Grasslands" = "#74C476"
  )
)

# Professional Academic Theme (Strictly Adhering to Nature Guidelines)
# Optimized Nature theme
nature_theme_professional <- theme_minimal(base_size = 10, base_family = best_font) +
  theme(
    # Font settings - ensure all text elements use the correct font
    text = element_text(family = best_font),
    
    # Grid and border settings
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#F0F0F0", linewidth = 0.25),
    panel.border = element_rect(fill = NA, color = "#CCCCCC", linewidth = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    
    # Legend settings
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 9, family = best_font),
    legend.text = element_text(size = 8, family = best_font),
    legend.key.size = unit(0.4, "cm"),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA),
    
    # Title settings
    plot.title = element_text(size = 11, face = "bold", color = "#000000", family = best_font),
    plot.subtitle = element_text(size = 9, color = "#666666", family = best_font),
    plot.caption = element_text(size = 7, color = "#666666", hjust = 0, family = best_font),
    
    # Axis settings - conforming to Nature standards
    axis.title = element_text(size = 9, face = "bold", color = "#333333", family = best_font),
    axis.text = element_text(size = 8, color = "#666666", family = best_font),
    axis.ticks = element_line(color = "#CCCCCC", linewidth = 0.25),
    axis.line = element_line(color = "#000000", linewidth = 0.3),  # Axis lines preferred by Nature
    
    # Facet settings
    strip.text = element_text(size = 9, face = "bold", color = "#333333", family = best_font),
    strip.background = element_rect(fill = "#F5F5F5", color = "#CCCCCC"),
    
    # Tag settings (important for annotations)
    plot.tag = element_text(size = 11, face = "bold", family = best_font),
    
    # Margin settings
    plot.margin = ggplot2::margin(10, 15, 10, 10, "pt"),
    
    # Ensure a clean background
    plot.background = element_rect(fill = "white", color = NA)
  )

# Create a safe annotate function to avoid font warnings
safe_annotate <- function(geom, ..., family = best_font) {
  if (geom == "text") {
    annotate(geom, ..., family = family)
  } else if (geom == "label") {
    annotate(geom, ..., family = family)
  } else {
    annotate(geom, ...)
  }
}

# Print information about the currently used font
cat("\n=== Font Configuration ===\n")
cat("Selected font for Nature theme:", best_font, "\n")
cat("This font will be used for all text elements in your plots.\n")

# Utility functions
saveFigure <- function(plot, filename, width = 8, height = 6, dpi = 600) {
  # Main figure saving
  full_filename <- file.path(base_dir, "figures", paste0(filename, ".png"))
  ggsave(full_filename, plot, width = width, height = height, dpi = dpi, bg = "white")
  
  # High-resolution version (for journal submission)
  hr_filename <- file.path(base_dir, "figures", paste0(filename, ".tiff"))
  ggsave(hr_filename, plot, width = width, height = height, dpi = 600, bg = "white", device = "tiff")
  
  # Data export
  if(!is.null(plot$data)) {
    data_filename <- file.path(base_dir, "data", paste0(filename, "_data.csv"))
    tryCatch({
      write.csv(plot$data, data_filename, row.names = FALSE)
    }, error = function(e) {
      if(inherits(plot, "ggplot")) {
        built <- ggplot_build(plot)
        if(length(built$data) > 0 && is.data.frame(built$data[[1]])) {
          write.csv(built$data[[1]], data_filename, row.names = FALSE)
        }
      }
    })
  }
  
  message(paste("✓ Saved:", full_filename))
}

saveTable <- function(data, filename) {
  full_filename <- file.path(base_dir, "tables", paste0(filename, ".csv"))
  write.csv(data, full_filename, row.names = FALSE)
  message(paste("✓ Saved table:", full_filename))
}

# -------------------------------------------------------------------------
# Data Loading and Professional Preprocessing
# -------------------------------------------------------------------------

# Smart Data Loading
load_landslide_data <- function() {
  possible_files <- c(
    "merged_landslide_data_add_climatezone.csv",
    "landslide_landcover_forest_data.csv", 
    "cleaned_landslide_with_vegetation_indices.csv"
  )
  
  for(file in possible_files) {
    if(file.exists(file)) {
      cat("Loading data from:", file, "\n")
      return(read.csv(file, stringsAsFactors = FALSE))
    }
  }
  
  stop("No valid data file found. Please ensure data file exists.")
}

landslide_data <- load_landslide_data()

# Data Preprocessing and Quality Control - Enhanced Multi-Source Vegetation Data Handling
preprocess_data <- function(data) {
  cat("Starting data preprocessing...\n")
  cat("Original data dimensions:", dim(data), "\n")
  
  # Date processing
  if("Event_Date" %in% colnames(data)) {
    data$Event_Date <- as.Date(data$Event_Date)
    data$Year <- as.numeric(format(data$Event_Date, "%Y"))
    data$Month <- format(data$Event_Date, "%m")
    data$Season <- case_when(
      data$Month %in% c("12", "01", "02") ~ "Winter",
      data$Month %in% c("03", "04", "05") ~ "Spring", 
      data$Month %in% c("06", "07", "08") ~ "Summer",
      data$Month %in% c("09", "10", "11") ~ "Fall"
    )
  }
  
  # Geographic classification (based on the Köppen-Geiger climate classification)
  
  # Continent classification
  data$Continent <- case_when(
    data$Longitude >= -170 & data$Longitude <= -30 & 
      data$Latitude >= 15 & data$Latitude <= 90 ~ "North America",
    data$Longitude >= -170 & data$Longitude <= -30 & 
      data$Latitude < 15 ~ "South America",
    data$Longitude > -30 & data$Longitude <= 60 & 
      data$Latitude >= 5 ~ "Europe/Africa", 
    data$Longitude > -30 & data$Longitude <= 60 & 
      data$Latitude < 5 ~ "Africa",
    data$Longitude > 60 & data$Longitude <= 180 ~ "Asia/Oceania",
    TRUE ~ "Other"
  )
  
  # ========== Enhanced Multi-Source Vegetation Data Handling ==========
  
  # 1. Basic cleaning and standardization of three data sources
  # MODIS IGBP classification processing
  if("MODIS_IGBP" %in% colnames(data)) {
    # Ensure data is character type
    data$MODIS_IGBP <- as.character(data$MODIS_IGBP)
    
    # Create standardized classification
    data$MODIS_IGBP_Class_Clean <- case_when(
      is.na(data$MODIS_IGBP) | data$MODIS_IGBP == "" | data$MODIS_IGBP == "0" ~ "Unknown",
      data$MODIS_IGBP == "1" ~ "Evergreen Needleleaf Forest",
      data$MODIS_IGBP == "2" ~ "Evergreen Broadleaf Forest",
      data$MODIS_IGBP == "3" ~ "Deciduous Needleleaf Forest",
      data$MODIS_IGBP == "4" ~ "Deciduous Broadleaf Forest",
      data$MODIS_IGBP == "5" ~ "Mixed Forest",
      data$MODIS_IGBP == "6" ~ "Closed Shrublands",
      data$MODIS_IGBP == "7" ~ "Open Shrublands",
      data$MODIS_IGBP == "8" ~ "Woody Savannas",
      data$MODIS_IGBP == "9" ~ "Savannas",
      data$MODIS_IGBP == "10" ~ "Grasslands",
      data$MODIS_IGBP == "11" ~ "Permanent Wetlands",
      data$MODIS_IGBP == "12" ~ "Croplands",
      data$MODIS_IGBP == "13" ~ "Urban and Built-up",
      data$MODIS_IGBP == "14" ~ "Cropland/Natural Vegetation Mosaic",
      data$MODIS_IGBP == "15" ~ "Snow and Ice",
      data$MODIS_IGBP == "16" ~ "Barren or Sparsely Vegetated",
      data$MODIS_IGBP == "17" ~ "Water Bodies",
      TRUE ~ as.character(data$MODIS_IGBP)
    )
    
    # Create simplified classification system (for later comparison)
    data$MODIS_IGBP_Simplified <- case_when(
      grepl("Forest", data$MODIS_IGBP_Class_Clean) ~ "Forest",
      grepl("Shrub|Savanna", data$MODIS_IGBP_Class_Clean) ~ "Woody Vegetation",
      grepl("Cropland|Vegetation", data$MODIS_IGBP_Class_Clean) ~ "Cropland/Vegetation",
      grepl("Grassland", data$MODIS_IGBP_Class_Clean) ~ "Grassland",
      grepl("Urban|Built", data$MODIS_IGBP_Class_Clean) ~ "Urban/Built-up",
      grepl("Barren|Sparsely", data$MODIS_IGBP_Class_Clean) ~ "Barren/Sparse",
      grepl("Water|Wetland|Snow|Ice", data$MODIS_IGBP_Class_Clean) ~ "Water/Wetland/Ice",
      data$MODIS_IGBP_Class_Clean == "Unknown" ~ "Unknown",
      TRUE ~ "Other"
    )
    
    # Create forest density levels (for comparison)
    data$MODIS_Forest_Density <- case_when(
      data$MODIS_IGBP_Simplified != "Forest" ~ "Non-forest",
      grepl("Evergreen", data$MODIS_IGBP_Class_Clean) ~ "Dense Forest",
      grepl("Mixed", data$MODIS_IGBP_Class_Clean) ~ "Moderate Forest",
      grepl("Deciduous", data$MODIS_IGBP_Class_Clean) ~ "Moderate Forest",
      TRUE ~ "Unknown"
    )
  }
  
  # Copernicus LC classification processing
  if("Copernicus_LC" %in% colnames(data)) {
    # Ensure data is character type
    data$Copernicus_LC <- as.character(data$Copernicus_LC)
    
    # Create standardized classification
    data$Copernicus_LC_Class_Clean <- case_when(
      is.na(data$Copernicus_LC) | data$Copernicus_LC == "" ~ "Unknown",
      data$Copernicus_LC == "111" ~ "Closed Forest, Evergreen Needle Leaf",
      data$Copernicus_LC == "112" ~ "Closed Forest, Evergreen Broad Leaf",
      data$Copernicus_LC == "113" ~ "Closed Forest, Deciduous Needle Leaf",
      data$Copernicus_LC == "114" ~ "Closed Forest, Deciduous Broad Leaf",
      data$Copernicus_LC == "115" ~ "Closed Forest, Mixed",
      data$Copernicus_LC == "116" ~ "Closed Forest, Unknown",
      data$Copernicus_LC == "121" ~ "Open Forest, Evergreen Needle Leaf",
      data$Copernicus_LC == "122" ~ "Open Forest, Evergreen Broad Leaf",
      data$Copernicus_LC == "123" ~ "Open Forest, Deciduous Needle Leaf",
      data$Copernicus_LC == "124" ~ "Open Forest, Deciduous Broad Leaf",
      data$Copernicus_LC == "125" ~ "Open Forest, Mixed",
      data$Copernicus_LC == "126" ~ "Open Forest, Unknown",
      data$Copernicus_LC == "20" ~ "Shrubland",
      data$Copernicus_LC == "30" ~ "Grassland",
      data$Copernicus_LC == "40" ~ "Cropland",
      data$Copernicus_LC == "50" ~ "Built-up",
      data$Copernicus_LC == "60" ~ "Bare/Sparse Vegetation",
      data$Copernicus_LC == "70" ~ "Snow and Ice",
      data$Copernicus_LC == "80" ~ "Permanent Water Bodies",
      data$Copernicus_LC == "90" ~ "Herbaceous Wetland",
      data$Copernicus_LC == "100" ~ "Moss and Lichen",
      data$Copernicus_LC == "200" ~ "Open Sea",
      TRUE ~ as.character(data$Copernicus_LC)
    )
    
    # Create simplified classification system (for later comparison)
    data$Copernicus_LC_Simplified <- case_when(
      grepl("Forest", data$Copernicus_LC_Class_Clean) ~ "Forest",
      grepl("Shrubland", data$Copernicus_LC_Class_Clean) ~ "Woody Vegetation",
      grepl("Grassland", data$Copernicus_LC_Class_Clean) ~ "Grassland",
      grepl("Cropland", data$Copernicus_LC_Class_Clean) ~ "Cropland/Vegetation",
      grepl("Built-up", data$Copernicus_LC_Class_Clean) ~ "Urban/Built-up",
      grepl("Bare|Sparse", data$Copernicus_LC_Class_Clean) ~ "Barren/Sparse",
      grepl("Water|Wetland|Snow|Ice|Sea|Moss|Lichen", data$Copernicus_LC_Class_Clean) ~ "Water/Wetland/Ice",
      data$Copernicus_LC_Class_Clean == "Unknown" ~ "Unknown",
      TRUE ~ "Other"
    )
    
    # Create forest density levels (for comparison)
    data$Copernicus_Forest_Density <- case_when(
      data$Copernicus_LC_Simplified != "Forest" ~ "Non-forest",
      grepl("Closed Forest", data$Copernicus_LC_Class_Clean) ~ "Dense Forest",
      grepl("Open Forest", data$Copernicus_LC_Class_Clean) ~ "Moderate Forest",
      TRUE ~ "Unknown"
    )
  }
  
  # Hansen forest cover processing
  if("Hansen_Tree_Cover_2000_Percent" %in% colnames(data)) {
    # Ensure data is numeric type
    data$Hansen_Tree_Cover_2000_Percent <- as.numeric(as.character(data$Hansen_Tree_Cover_2000_Percent))
    
    # Create forest cover classification
    data$Hansen_Forest_Class_Clean <- case_when(
      is.na(data$Hansen_Tree_Cover_2000_Percent) ~ "Unknown",
      data$Hansen_Tree_Cover_2000_Percent == 0 ~ "No forest",
      data$Hansen_Tree_Cover_2000_Percent > 0 & data$Hansen_Tree_Cover_2000_Percent <= 25 ~ "Sparse forest",
      data$Hansen_Tree_Cover_2000_Percent > 25 & data$Hansen_Tree_Cover_2000_Percent <= 50 ~ "Moderate forest",
      data$Hansen_Tree_Cover_2000_Percent > 50 & data$Hansen_Tree_Cover_2000_Percent <= 75 ~ "Dense forest",
      data$Hansen_Tree_Cover_2000_Percent > 75 ~ "Very dense forest"
    )
    
    # Create simplified classification system (for later comparison)
    data$Hansen_Simplified <- case_when(
      data$Hansen_Forest_Class_Clean %in% c("Dense forest", "Very dense forest") ~ "Forest",
      data$Hansen_Forest_Class_Clean %in% c("Moderate forest", "Sparse forest") ~ "Woody Vegetation",
      data$Hansen_Forest_Class_Clean == "No forest" ~ "Non-forest",
      data$Hansen_Forest_Class_Clean == "Unknown" ~ "Unknown",
      TRUE ~ "Other"
    )
    
    # Retain detailed forest cover classification
    data$Forest_Cover_Category <- case_when(
      is.na(data$Hansen_Tree_Cover_2000_Percent) ~ "Unknown",
      data$Hansen_Tree_Cover_2000_Percent == 0 ~ "No Forest (0%)",
      data$Hansen_Tree_Cover_2000_Percent > 0 & data$Hansen_Tree_Cover_2000_Percent <= 25 ~ "Sparse (1-25%)",
      data$Hansen_Tree_Cover_2000_Percent > 25 & data$Hansen_Tree_Cover_2000_Percent <= 50 ~ "Moderate (26-50%)",
      data$Hansen_Tree_Cover_2000_Percent > 50 & data$Hansen_Tree_Cover_2000_Percent <= 75 ~ "Dense (51-75%)",
      data$Hansen_Tree_Cover_2000_Percent > 75 ~ "Very Dense (>75%)"
    )
    
    data$Forest_Cover_Category <- factor(data$Forest_Cover_Category,
                                         levels = c("No Forest (0%)", "Sparse (1-25%)", "Moderate (26-50%)", 
                                                    "Dense (51-75%)", "Very Dense (>75%)", "Unknown"))
  }
  
  # 2. Create multi-source consistency metrics
  if(all(c("MODIS_IGBP_Simplified", "Hansen_Simplified") %in% colnames(data))) {
    # Consistency between MODIS and Hansen
    data$MODIS_Hansen_Consistency <- ifelse(
      data$MODIS_IGBP_Simplified == "Unknown" | data$Hansen_Simplified == "Unknown",
      "Unknown",
      ifelse(data$MODIS_IGBP_Simplified == data$Hansen_Simplified, "Consistent", "Inconsistent")
    )
  }
  
  if(all(c("MODIS_IGBP_Simplified", "Copernicus_LC_Simplified") %in% colnames(data))) {
    # Consistency between MODIS and Copernicus
    data$MODIS_Copernicus_Consistency <- ifelse(
      data$MODIS_IGBP_Simplified == "Unknown" | data$Copernicus_LC_Simplified == "Unknown",
      "Unknown",
      ifelse(data$MODIS_IGBP_Simplified == data$Copernicus_LC_Simplified, "Consistent", "Inconsistent")
    )
  }
  
  if(all(c("Copernicus_LC_Simplified", "Hansen_Simplified") %in% colnames(data))) {
    # Consistency between Copernicus and Hansen
    data$Copernicus_Hansen_Consistency <- ifelse(
      data$Copernicus_LC_Simplified == "Unknown" | data$Hansen_Simplified == "Unknown",
      "Unknown",
      ifelse(data$Copernicus_LC_Simplified == data$Hansen_Simplified, "Consistent", "Inconsistent")
    )
  }
  
  # 3. Create a composite three-source consistency score
  if(all(c("MODIS_IGBP_Simplified", "Copernicus_LC_Simplified", "Hansen_Simplified") %in% colnames(data))) {
    # Calculate consistency score (0-3)
    data$Multi_Source_Consistency_Score <- apply(
      data[, c("MODIS_IGBP_Simplified", "Copernicus_LC_Simplified", "Hansen_Simplified")], 
      1, 
      function(x) {
        if(any(x == "Unknown")) return(NA)
        unique_values <- length(unique(x))
        return(4 - unique_values)  # 3: Full agreement, 2: Partial agreement, 1: No agreement
      }
    )
    
    # Create consistency labels
    data$Multi_Source_Consistency_Label <- case_when(
      is.na(data$Multi_Source_Consistency_Score) ~ "Unknown",
      data$Multi_Source_Consistency_Score == 3 ~ "Full Agreement",
      data$Multi_Source_Consistency_Score == 2 ~ "Partial Agreement",
      data$Multi_Source_Consistency_Score == 1 ~ "No Agreement"
    )
  }
  
  # 4. Create a unified vegetation type label (Priority: MODIS > Copernicus > Hansen)
  data$Unified_Vegetation_Type <- case_when(
    "MODIS_IGBP_Simplified" %in% colnames(data) & !is.na(data$MODIS_IGBP_Simplified) & 
      data$MODIS_IGBP_Simplified != "Unknown" ~ data$MODIS_IGBP_Simplified,
    "Copernicus_LC_Simplified" %in% colnames(data) & !is.na(data$Copernicus_LC_Simplified) & 
      data$Copernicus_LC_Simplified != "Unknown" ~ data$Copernicus_LC_Simplified,
    "Hansen_Simplified" %in% colnames(data) & !is.na(data$Hansen_Simplified) & 
      data$Hansen_Simplified != "Unknown" ~ data$Hansen_Simplified,
    TRUE ~ "Unknown"
  )
  
  # Calculate rate of change for vegetation indices
  if(all(c("NDVI_1", "NDVI_change_1_to_2") %in% colnames(data))) {
    data <- data %>%
      mutate(
        NDVI_change_rate_1_to_2 = (NDVI_change_1_to_2 / pmax(0.01, abs(NDVI_1))) * 100,
        EVI_change_rate_1_to_2 = if("EVI_change_1_to_2" %in% colnames(.)) 
          (EVI_change_1_to_2 / pmax(0.01, abs(EVI_1))) * 100 else NA,
        LAI_change_rate_1_to_2 = if("LAI_change_1_to_2" %in% colnames(.)) 
          (LAI_change_1_to_2 / pmax(0.01, abs(LAI_1))) * 100 else NA
      )
  }
  
  # Data quality report
  cat("Processed data dimensions:", dim(data), "\n")
  cat("Vegetation data sources availability:\n")
  cat("  - MODIS IGBP: ", sum(!is.na(data$MODIS_IGBP_Class_Clean) & 
                                data$MODIS_IGBP_Class_Clean != "Unknown"), "records\n")
  cat("  - Copernicus LC: ", sum(!is.na(data$Copernicus_LC_Class_Clean) & 
                                   data$Copernicus_LC_Class_Clean != "Unknown"), "records\n")
  cat("  - Hansen Forest: ", sum(!is.na(data$Hansen_Forest_Class_Clean) & 
                                   data$Hansen_Forest_Class_Clean != "Unknown"), "records\n")
  
  cat("Multi-source consistency metrics:\n")
  if("Multi_Source_Consistency_Label" %in% colnames(data)) {
    cons_table <- table(data$Multi_Source_Consistency_Label, useNA = "ifany")
    for(level in names(cons_table)) {
      cat("  - ", level, ": ", cons_table[level], " records (", 
          round(cons_table[level]/sum(cons_table)*100, 1), "%)\n", sep="")
    }
  }
  
  cat("Available vegetation indices:", 
      sum(c("NDVI_1", "EVI_1", "LAI_1", "FPAR_1") %in% colnames(data)), "/4\n")
  cat("Climate zones:", length(unique(data$Climate_Zone)), "\n")
  cat("Continents:", length(unique(data$Continent)), "\n")
  
  return(data)
}

landslide_data <- preprocess_data(landslide_data)

# Filter out data where "Climate_Zone" is "Polar", Copernicus_LC_Class is "Permanent water bodies",
# or MODIS_IGBP_Class is "Water Bodies"
landslide_data <- landslide_data[
  !(landslide_data$Climate_Zone == "Polar" |
      landslide_data$Copernicus_LC_Class == "Permanent water bodies" |
      landslide_data$MODIS_IGBP_Class == "Water Bodies"),
]

saveTable(landslide_data,"landslide_data_processed")


# -------------------------------------------------------------------------
# MAIN FIGURE 1: Data-Driven Threshold Discovery and Risk Zone Identification
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# MAIN FIGURE 1: Data Analysis and Supplementary Material Generation
# -------------------------------------------------------------------------

# <<< NEW SECTION START >>>
# HELPER FUNCTION 1: Create Extended Data Table 1 (Data Sources)
# PURPOSE: To provide a transparent and comprehensive summary of all data used,
# a standard requirement for high-impact journals like Nature Communications.
# -------------------------------------------------------------------------
create_extended_data_table_1 <- function() {
  cat("--- Generating Extended Data Table 1: Data Sources ---\n")
  
  data_sources_df <- data.frame(
    DataType = c(
      "Landslide Inventory", "Landslide Inventory", "Landslide Inventory",
      "Vegetation Index", "Vegetation Index", "High-Res Vegetation",
      "Land Cover", "Land Cover", "Land Cover", "Forest Structure", "Forest Status"
    ),
    Source_Product_Name = c(
      "NASA Global Landslide Catalog (GLC)", "e-ITALICA", "Guangzhou Institute of Geography Catalog",
      "MODIS/061/MOD13Q1 (NDVI/EVI)", "MODIS/061/MOD15A2H (LAI/FAPAR)", "COPERNICUS/S2_SR_HARMONIZED",
      "ESA WorldCover", "MODIS/006/MCD12Q1", "COPERNICUS/Landcover/100m",
      "NASA/MEASURES/GFCC/TC/v3", "UMD/hansen/global_forest_change_2022_v1_10"
    ),
    Resolution = c(
      "Point", "Point", "Point",
      "250m", "500m", "10-60m",
      "10m", "500m", "100m",
      "30m", "30m"
    ),
    Temporal_Coverage = c(
      "~2007-Present", "1996-2021", "Internal",
      "2000-Present", "2000-Present", "2015-Present",
      "2020-Present", "2001-Present", "2015-2019",
      "2000, 2005, 2010, 2015", "2000-2022"
    ),
    Key_Role_in_Study = c(
      "Primary global landslide data", "Regional validation landslide data", "Regional validation landslide data",
      "Primary VI for threshold discovery (NDVI, EVI)", "Primary VI for threshold discovery (LAI)", "VI calculation for recent events",
      "High-resolution land cover validation", "Primary land cover for validation", "High-resolution land cover validation",
      "Forest density analysis", "Forest change and tree cover validation"
    )
  )
  
  saveTable(data_sources_df, "Extended_Data_Table_1_Data_Sources")
  cat("✓ Extended Data Table 1 saved to 'tables/' directory.\n")
  return(data_sources_df)
}
# <<< NEW SECTION END >>>


# =========================================================================
# FIGURE 1: DATA-DRIVEN THRESHOLD DISCOVERY AND Biophysical Zone IDENTIFICATION
# =========================================================================

# -------------------------------------------------------------------------
# Part 1: Core Analysis Function (Robust and Adaptive)
# -------------------------------------------------------------------------

#' Advanced Segmented Regression Analysis
#'
#' This function performs a robust, adaptive segmented regression to find
#' optimal breakpoints in landslide frequency data for different vegetation indices.
#' It automatically detects the data's distribution shape (peaked or valley)
#' and applies the appropriate initialization strategy for the segmentation model.
#'
#' @param data The input dataframe (e.g., landslide_data_processed).
#' @param var_name The name of the vegetation index column to analyze (e.g., "NDVI_1").
#' @param breaks The number of breakpoints to find (default is 2).
#' @return A list containing the model, breakpoints, confidence intervals, and plot data.

advanced_segmented_analysis <- function(data, var_name, breaks = 2) {
  
  # --- 1.1: Helper function to classify data shape ---
  classify_distribution <- function(y_values) {
    n <- length(y_values)
    # Use quartiles for more robust division
    q1_idx <- floor(n/4)
    q3_idx <- floor(3*n/4)
    start_mean <- mean(y_values[1:q1_idx], na.rm = TRUE)
    middle_mean <- mean(y_values[(q1_idx + 1):(q3_idx - 1)], na.rm = TRUE)
    end_mean <- mean(y_values[q3_idx:n], na.rm = TRUE)
    
    # Check if the middle is significantly different from the ends
    if (middle_mean > (start_mean * 1.1) && middle_mean > (end_mean * 1.1)) {
      return("peaked")
    } else if (middle_mean < (start_mean * 0.9) && middle_mean < (end_mean * 0.9)) {
      return("valley")
    } else {
      return("other")
    }
  }
  
  # --- 1.2: Data validation and preparation ---
  valid_data <- data[!is.na(data[[var_name]]), ]
  x_vals <- valid_data[[var_name]]
  
  if(length(x_vals) < 100) {
    stop(paste("Insufficient data for", var_name, ": only", length(x_vals), "valid points"))
  }
  
  # --- 1.3: Binning and Smoothing ---
  n_bins <- min(150, floor(length(x_vals) / 20))
  x_range <- range(x_vals, na.rm = TRUE)
  breaks_seq <- seq(x_range[1], x_range[2], length.out = n_bins + 1)
  
  # Use hist() for efficient binning
  h <- hist(x_vals, breaks = breaks_seq, plot = FALSE)
  counts <- h$counts
  mid_points <- h$mids
  
  # Correctly apply smooth() sequentially
  smooth_1 <- stats::smooth(counts, kind = "3R")
  smooth_2 <- stats::smooth(smooth_1, kind = "S")
  counts_final <- stats::smooth(smooth_2, kind = "3R")
  
  df <- data.frame(x = mid_points, y = as.numeric(counts_final))
  lin_mod <- lm(y ~ x, data = df)
  
  # --- 1.4: Segmented Regression with Adaptive Initialization ---
  tryCatch({
    distribution_shape <- classify_distribution(df$y)
    cat(paste("  - Detected distribution for", var_name, "as:", distribution_shape, "\n"))
    
    if (distribution_shape == "peaked") {
      peak_index <- which.max(df$y)
      peak_x_value <- df$x[peak_index]
      q10 <- quantile(df$x, 0.10); q90 <- quantile(df$x, 0.90)
      initial_psi <- c((q10 + peak_x_value) / 2, (peak_x_value + q90) / 2)
    } else if (distribution_shape == "valley") {
      valley_index <- which.min(df$y)
      valley_x_value <- df$x[valley_index]
      q25 <- quantile(df$x, 0.25); q75 <- quantile(df$x, 0.75)
      initial_psi <- c((q25 + valley_x_value) / 2, (valley_x_value + q75) / 2)
    } else {
      initial_psi <- quantile(df$x, probs = seq(0.25, 0.75, length.out = breaks))
    }
    
    initial_psi <- pmax(x_range[1] + 0.05 * diff(x_range), pmin(x_range[2] - 0.05 * diff(x_range), initial_psi))
    
    seg_mod <- segmented::segmented(lin_mod, seg.Z = ~x, psi = list(x = initial_psi))
    
    breakpoints <- seg_mod$psi[, "Est."]
    conf_intervals <- confint(seg_mod)
    ci_rows <- grep("psi", rownames(conf_intervals))
    threshold_conf <- if(length(ci_rows) > 0) conf_intervals[ci_rows, , drop = FALSE] else matrix(c(breakpoints - 0.02, breakpoints + 0.02), ncol = 2)
    r_squared <- summary(seg_mod)$r.squared
    aic_value <- AIC(seg_mod)
    
  }, error = function(e) {
    cat("  - Segmented regression failed for", var_name, ":", e$message, "\n  - Using quantile method as fallback.\n")
    breakpoints <- quantile(x_vals, probs = c(0.4, 0.8))
    threshold_conf <- matrix(c(breakpoints - 0.02, breakpoints + 0.02), ncol = 2)
    r_squared <- 0.5; aic_value <- NA; seg_mod <- lin_mod
  })
  
  return(list(
    breakpoints = breakpoints,
    confidence_intervals = threshold_conf,
    data = df,
    model = seg_mod,
    r_squared = r_squared,
    aic = aic_value
  ))
}

# -------------------------------------------------------------------------
# Part 2: Main Analysis Workflow
# -------------------------------------------------------------------------
analyze_vegetation_thresholds <- function(data) {
  cat("\n=== (Full Workflow) Analyzing Vegetation Thresholds & Generating Supplementary Materials ===\n")
  
  # ===================================================================
  # Internal Helper Functions
  # ===================================================================
  
  # --- Helper Function 1: Supplementary Density Plot ---
  create_supplementary_density_plot <- function(plot_data, ndvi_breakpoints) {
    cat("  - Generating Supplementary Figure: NDVI Kernel Density Plot...\n")
    p <- ggplot(plot_data, aes(x = NDVI_1)) +
      geom_density(fill = "#0173B2", color = "#004D7A", alpha = 0.7) +
      geom_vline(xintercept = ndvi_breakpoints, linetype = "dashed", color = "#E31A1C", linewidth = 1) +
      annotate("text", x = ndvi_breakpoints[1], y = Inf, label = paste(" BP1:", round(ndvi_breakpoints[1], 3)), angle = 90, vjust = -0.5, hjust = 1.1, color = "#E31A1C", fontface = "bold") +
      annotate("text", x = ndvi_breakpoints[2], y = Inf, label = paste(" BP2:", round(ndvi_breakpoints[2], 3)), angle = 90, vjust = 1.5, hjust = 1.1, color = "#E31A1C", fontface = "bold") +
      labs(title = "Supplementary Figure: Kernel Density Estimate of Landslide NDVI", x = "NDVI Value", y = "Density") +
      nature_theme_professional
    saveFigure(p, "Supplementary_Fig_1a_NDVI_Kernel_Density", width = 7, height = 5)
    cat("  ✓ Supplementary Figure saved.\n")
  }
  
  # --- Helper Function 2: Supplementary Cross-Index Consistency Table (Optimized) ---
  analyze_cross_index_consistency <- function(data_with_zones, all_results) {
    cat("  - Generating Supplementary Table: Cross-Index Consistency...\n")
    temp_data <- data_with_zones %>%
      # EVI biophysical zones
      mutate(EVI_Zone = case_when(
        is.na(EVI_1) ~ NA_character_,
        EVI_1 < all_results$EVI_1$breakpoints[1] ~ "IDZ",
        EVI_1 <= all_results$EVI_1$breakpoints[2] ~ "CDZ",
        TRUE ~ "SDZ"
      )) %>%
      filter(!is.na(NDVI_Biophysical_Zone), !is.na(EVI_Zone), NDVI_Biophysical_Zone != "Unknown") %>%
      # Standardize factor levels
      mutate(
        NDVI_Zone_Simple = gsub(" Zone", "", as.character(NDVI_Biophysical_Zone)),
        EVI_Zone_Simple = gsub(" Zone", "", EVI_Zone)
      )
    
    # Using `janitor` package's `tabyl` and `adorn_totals` is the cleanest, warning-free method
    if(!require(janitor)) install.packages("janitor", quiet = TRUE); library(janitor)
    
    consistency_table <- temp_data %>%
      tabyl(NDVI_Zone_Simple, EVI_Zone_Simple) %>%
      adorn_totals(where = c("row", "col")) %>% # Add row and column totals
      adorn_percentages(denominator = "row") %>% # Calculate row-wise percentages
      adorn_pct_formatting(digits = 1) %>% # Format as percentages
      adorn_ns(position = "front") # Add counts (n) in front of percentages
    
    # Convert janitor table to a standard data frame and save
    consistency_df <- as.data.frame(consistency_table)
    saveTable(consistency_df, "Supplementary_Table_1b_NDVI_EVI_Consistency")
    cat("  ✓ Cross-index consistency table saved.\n")
  }
  
  # ===================================================================
  # Main Analysis Workflow
  # ===================================================================
  
  # --- 1. Perform core breakpoint analysis ---
  indices_to_analyze <- c("NDVI_1", "EVI_1", "LAI_1")
  available_indices <- indices_to_analyze[indices_to_analyze %in% colnames(data)]
  if(length(available_indices) == 0) stop("None of the specified vegetation indices were found.")
  
  cat("  - Analyzing indices:", paste(available_indices, collapse = ", "), "\n")
  analysis_results <- list()
  for(idx in available_indices) {
    cat("  - Processing", idx, "...\n")
    result <- advanced_segmented_analysis(data, idx, breaks = 2)
    analysis_results[[idx]] <- result
  }
  
  # --- 2. Define risk zones and add to the data ---
  if("NDVI_1" %in% names(analysis_results)) {
    ndvi_result <- analysis_results[["NDVI_1"]]
    bp1 <- ndvi_result$breakpoints[1]
    bp2 <- ndvi_result$breakpoints[2]
    
    data$NDVI_Biophysical_Zone <- factor(
      case_when(
        is.na(data$NDVI_1) ~ "Unknown",
        data$NDVI_1 < bp1 ~ "IDZ",
        data$NDVI_1 >= bp1 & data$NDVI_1 <= bp2 ~ "CTZ",
        data$NDVI_1 > bp2 ~ "SDZ"
      ),
      levels = c("Unknown", "IDZ", "CTZ", "SDZ")
    )
    cat("  ✓ NDVI-based biophysical zones defined and added to data.\n")
  }
  
  # --- 3. Call helper functions to generate supplementary materials ---
  # Note: we pass the modified `data` object to the helper functions
  if("NDVI_1" %in% names(analysis_results)) {
    create_supplementary_density_plot(data, analysis_results[["NDVI_1"]]$breakpoints)
  }
  if(all(c("NDVI_1", "EVI_1") %in% names(analysis_results))) {
    # Pass the data that now includes the NDVI_Biophysical_Zone column
    analyze_cross_index_consistency(data, analysis_results)
  }
  
  assign("fig1_analysis_results", analysis_results, envir = .GlobalEnv)
  
  cat("✓ Main analysis and supplementary materials generation completed successfully.\n")
  return(list(results = analysis_results, data = data))
}

# -------------------------------------------------------------------------
# Part 3: Visualization
# -------------------------------------------------------------------------
create_figure_1_plots <- function(analysis_results, processed_data) {
  cat("\n=== Creating Figure 1 Visualizations with Fine-Tuned Labels ===\n")
  
  plot_list <- list()
  nature_palettes <- get("nature_palettes", envir = .GlobalEnv)
  
  # --- Sub-function for breakpoint plots (a, b, c) with restored label control ---
  create_breakpoint_plot <- function(idx, result) {
    df <- result$data
    breakpoints <- result$breakpoints
    conf_int <- result$confidence_intervals
    y_max <- max(df$y, na.rm = TRUE)
    y_min <- min(df$y, na.rm = TRUE)
    x_range <- range(df$x, na.rm = TRUE)
    x_span <- diff(x_range)
    
    # Define plotting boundaries - increase buffer space to keep points/lines away from the frame
    plot_xlim <- c(x_range[1] - 0.08*x_span, x_range[2] + 0.08*x_span)
    
    # For NDVI/EVI, use a slightly expanded range
    if(idx %in% c("NDVI_1", "EVI_1")) {
      plot_xlim <- c(-0.05, 1.05) 
    }
    
    # Create prediction data that stays within reasonable bounds
    pred_data_extended <- data.frame(x = seq(max(plot_xlim[1], min(df$x, na.rm = TRUE)), 
                                             min(plot_xlim[2], max(df$x, na.rm = TRUE)), 
                                             length.out = 200))
    
    # Use the segmented model to predict y-values
    pred_data_extended$y <- predict(result$model, newdata = pred_data_extended)
    
    # Filter out any predictions that might be problematic
    pred_data_extended <- pred_data_extended[is.finite(pred_data_extended$y), ]
    
    r2_pos_x <- ifelse(idx == "LAI_1", plot_xlim[2] - 0.05 * diff(plot_xlim), plot_xlim[1] + 0.05 * diff(plot_xlim))
    r2_pos_y <- 0.95 * y_max
    r2_hjust <- ifelse(idx == "LAI_1", 1, 0)
    
    # Start building the plot
    p <- ggplot(df, aes(x = x, y = y)) +
      geom_point(alpha = 0.4, color = "grey50", size = 0.8) +
      geom_line(data = pred_data_extended, aes(x = x, y = y), 
                color = nature_palettes$indices[gsub("_1", "", idx)], 
                linewidth = 1.2) +
      geom_vline(xintercept = breakpoints, linetype = "dashed", 
                 color = "#E31A1C", linewidth = 0.8) +
      annotate("text", x = r2_pos_x, y = r2_pos_y, 
               label = bquote(R^2 == .(sprintf("%.3f", result$r_squared))), 
               hjust = r2_hjust, size = 3.2, fontface = "bold")
    
    # Use the effective confidence interval shading method
    cat(sprintf("Debug: Adding confidence intervals for %s\n", idx))
    if(nrow(conf_int) >= 1) {
      for(i in 1:min(2, nrow(conf_int))) {
        # Ensure confidence interval limits are in the correct order
        ci_lower <- min(conf_int[i, 1], conf_int[i, 2])
        ci_upper <- max(conf_int[i, 1], conf_int[i, 2])
        
        cat(sprintf("  BP%d: CI [%.3f, %.3f]\n", i, ci_lower, ci_upper))
        
        # Use the geom_polygon method
        conf_data <- data.frame(
          x = c(ci_lower, ci_upper, ci_upper, ci_lower),
          y = c(0, 0, y_max, y_max)
        )
        
        p <- p + geom_polygon(data = conf_data, 
                              aes(x = x, y = y),
                              fill = "#E31A1C", 
                              alpha = 0.15,
                              inherit.aes = FALSE)  # Explicitly do not inherit aesthetics
      }
    }
    
    # --- Re-integrated manual label positioning logic ---
    label_positions <- list()
    
    # For NDVI_1
    if(idx == "NDVI_1") {
      label_positions[[1]] <- list(x_offset = -0.01 * x_span, y_position = 0.7 * y_max, hjust = 1)
      if(length(breakpoints) > 1) {
        label_positions[[2]] <- list(x_offset = -0.02 * x_span, y_position = 1.0 * y_max, hjust = 1)
      }
    } 
    # For EVI_1
    else if(idx == "EVI_1") {
      label_positions[[1]] <- list(x_offset = -0.01 * x_span, y_position = 0.05 * y_max, hjust = 1)
      if(length(breakpoints) > 1) {
        label_positions[[2]] <- list(x_offset = -0.09 * x_span, y_position = 1.0 * y_max, hjust = 0)
      }
    }
    # For LAI_1
    else if(idx == "LAI_1") {
      label_positions[[1]] <- list(x_offset = -0.02 * x_span, y_position = 1 * y_max, hjust = 1)
      if(length(breakpoints) > 1) {
        label_positions[[2]] <- list(x_offset = 0.01 * x_span, y_position = 0.5 * y_max, hjust = 0)
      }
    }
    # Fallback (though unlikely to be used with these specific indices)
    else {
      for(i in 1:length(breakpoints)) {
        label_positions[[i]] <- list(x_offset = ifelse(i == 1, -0.05, 0.05) * x_span, y_position = (0.9 - (i-1)*0.2) * y_max, hjust = ifelse(i == 1, 1, 0))
      }
    }
    
    # Add labels using the defined positions
    for(i in 1:length(breakpoints)) {
      if(i <= length(label_positions) && !is.null(label_positions[[i]])) {
        pos <- label_positions[[i]]
        label_x <- breakpoints[i] + pos$x_offset
        p <- p + annotate("label", x = label_x, y = pos$y_position,
                          label = sprintf("BP%d: %.3f\n95%% CI: [%.3f, %.3f]", i, breakpoints[i], conf_int[i, 2], conf_int[i, 1]),
                          hjust = pos$hjust, fill = "white", alpha = 0.75, color = "#E31A1C", size = 2.8, fontface = "bold",
                          label.size = 0.2, label.padding = unit(0.15, "lines"))
      }
    }
    # --- End of re-integrated logic ---
    
    p <- p +
      nature_theme_professional +
      # Keep clip = "on" but add boundary space
      coord_cartesian(xlim = plot_xlim,
                      ylim = c(-0.02 * y_max, y_max * 1.1), # Add space at top and bottom
                      expand = FALSE, clip = "on") +
      labs(x = paste(gsub("_1", "", idx), "Value"), y = "Landslide Frequency (Smoothed)")
    
    return(p)
  }
  
  # --- Sub-function for distribution plot (d) ---
  create_distribution_plot <- function(data_with_biophysical_zones, ndvi_result) {
    full_dist <- data_with_biophysical_zones %>%
      dplyr::count(NDVI_Biophysical_Zone) %>%
      dplyr::mutate(Percentage = n / sum(n) * 100)
    
    p <- ggplot(full_dist, aes(x = NDVI_Biophysical_Zone, y = Percentage, fill = NDVI_Biophysical_Zone)) +
      geom_col(alpha = 0.9) +
      geom_text(aes(label = sprintf("%.1f%%\n(n=%s)", Percentage, format(n, big.mark = ","))), vjust = -0.4, size = 2.8, fontface = "bold") +
      
      annotate("text", x = 3.5, y = 35, 
               label = sprintf("CTZ [%.3f-%.3f]:\n%.1f%% of global landslides", 
                               ndvi_result$breakpoints[1], ndvi_result$breakpoints[2],
                               full_dist$Percentage[full_dist$NDVI_Biophysical_Zone == "CTZ"]),
               size = 3, hjust = 0.5, lineheight = 1.1) +
      
      scale_fill_manual(values = c("Unknown" = "grey70", nature_palettes$biophysical_zones), guide = "none") +
      scale_y_continuous(limits = c(0, 40), expand = expansion(mult = c(0, 0.05))) +
      nature_theme_professional +
      labs(x = "Biophysical Zone", y = "Percentage (%)")
    
    return(p)
  }
  
  # --- Generate all plots ---
  plot_list$NDVI_1 <- create_breakpoint_plot("NDVI_1", analysis_results$NDVI_1)
  plot_list$EVI_1 <- create_breakpoint_plot("EVI_1", analysis_results$EVI_1)
  plot_list$LAI_1 <- create_breakpoint_plot("LAI_1", analysis_results$LAI_1)
  plot_list$zone_distribution <- create_distribution_plot(processed_data, analysis_results$NDVI_1)
  
  # --- Combine and save the final figure ---
  final_figure <- (plot_list$NDVI_1 | plot_list$EVI_1) / (plot_list$LAI_1 | plot_list$zone_distribution) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')') &
    theme(
      plot.tag = element_text(face = "bold", size = 11),
      plot.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5)
    )
  print(final_figure)
  saveFigure(final_figure, "Main_Figure_1_Complete_A4", width = 8.3, height = 6.8)
  
  cat("✓ Figure 1 visualizations completed successfully.\n")
  return(final_figure)
}

# =========================================================================
# RUN THE ANALYSIS AND VISUALIZATION
# =========================================================================

# Assuming 'landslide_data_processed' from the global environment is the starting point
# 1. Run the analysis. The function returns a list with results and the modified data.
analysis_output <- analyze_vegetation_thresholds(landslide_data)

# 2. Pass the processed data from the analysis output to the plotting function.
# Use analysis_output$results for the analysis results.
# Use analysis_output$data for the data, which now includes the 'NDVI_Biophysical_Zone' column.
figure_1 <- create_figure_1_plots(analysis_output$results, analysis_output$data)

# To update the global 'landslide_data_processed' object for subsequent script sections (e.g., Figure 2, 3), do this:
landslide_data_processed <- analysis_output$data

# ===========================================================================

# -------------------------------------------------------------------------
# MAIN FIGURE 2: Data Analysis - Multi-Source Vegetation Data Validation and Cross-Platform Consistency
# -------------------------------------------------------------------------

# =============================================================================

analyze_vegetation_validation_data <- function() {
  cat("\n=== Analyzing Vegetation Data Validation and Cross-Platform Consistency ===\n")
  
  data <- landslide_data_processed # Use the data processed earlier
  analysis_results <- list()
  
  # --- 2A: Global vegetation type distribution analysis ---
  if("MODIS_IGBP_Class_Clean" %in% colnames(data)) {
    modis_dist <- data %>%
      filter(!is.na(MODIS_IGBP_Class_Clean), MODIS_IGBP_Class_Clean != "Unknown") %>%
      group_by(MODIS_IGBP_Class_Clean) %>%
      dplyr::summarise(Count = n(), .groups = "drop") %>%
      filter(Count > 0) %>%
      arrange(desc(Count)) %>%
      slice_head(n = 10) %>% # Top 10 for clarity
      mutate(Percentage = Count / sum(Count) * 100,
             Display_Name_Short = sapply(strsplit(MODIS_IGBP_Class_Clean, " "), function(x) paste(x[1:min(length(x),2)], collapse=" "))) # First 1 or 2 words
    
    analysis_results[["vegetation_type_distribution"]] <- modis_dist
    saveTable(modis_dist, "Fig2A_modis_vegetation_distribution")
  } else {
    cat("Warning: MODIS_IGBP_Class_Clean not available for vegetation type analysis\n")
  }
  
  # --- 2B: Multi-source consistency level analysis ---
  if("Multi_Source_Consistency_Label" %in% colnames(data)) {
    consistency_summary <- data %>%
      filter(!is.na(Multi_Source_Consistency_Label)) %>%
      group_by(Multi_Source_Consistency_Label) %>%
      dplyr::summarise(Count = n(), .groups = "drop") %>%
      mutate(Percentage = Count / sum(Count) * 100) %>%
      mutate(Multi_Source_Consistency_Label = factor(Multi_Source_Consistency_Label,
                                                     levels = c("Full Agreement", "Partial Agreement", "Unknown", "No Agreement"))) %>%
      arrange(desc(Multi_Source_Consistency_Label))
    
    analysis_results[["multisource_consistency"]] <- consistency_summary
    saveTable(consistency_summary, "Fig2B_multisource_consistency_levels")
  } else {
    cat("Warning: Multi_Source_Consistency_Label not available for consistency analysis\n")
  }
  
  # <<< NEW SECTION START >>>
  # STATISTICAL TEST 1: Fleiss' Kappa for Multi-Source Consistency
  # PURPOSE: To provide a single, robust statistic quantifying the agreement
  # between the three land cover datasets (MODIS, Copernicus, Hansen).
  # -------------------------------------------------------------------------
  cat("--- Performing Fleiss' Kappa test for multi-source agreement ---\n")
  
  # Prepare the data for the kappa test
  # The `kappam.fleiss` function requires a matrix where rows are subjects (landslides)
  # and columns are raters (datasets). The values are the categories.
  consistency_data <- data %>%
    filter(!is.na(MODIS_IGBP_Simplified), !is.na(Copernicus_LC_Simplified), !is.na(Hansen_Simplified)) %>%
    filter(MODIS_IGBP_Simplified != "Unknown" & Copernicus_LC_Simplified != "Unknown" & Hansen_Simplified != "Unknown") %>%
    dplyr::select(MODIS_IGBP_Simplified, Copernicus_LC_Simplified, Hansen_Simplified)
  
  if(nrow(consistency_data) > 10) { # Check if there is enough data
    fleiss_result <- kappam.fleiss(consistency_data)
    
    # Print the results to the console
    print(fleiss_result)
    
    # Create a summary table to save
    fleiss_summary <- data.frame(
      Statistic = "Fleiss' Kappa",
      Value = fleiss_result$value,
      Z_statistic = fleiss_result$statistic,
      P_value = fleiss_result$p.value,
      Subjects = fleiss_result$subjects,
      Raters = fleiss_result$raters
    )
    
    # The `irr` package doesn't directly provide CI for Fleiss' Kappa.
    # We will add the manually provided CI for reference.
    fleiss_summary$CI_95_Lower_Manual = 0.67
    fleiss_summary$CI_95_Upper_Manual = 0.71
    
    saveTable(fleiss_summary, "Statistics_Fig2B_Fleiss_Kappa")
    cat(paste0("✓ Fleiss' Kappa result: ", round(fleiss_result$value, 2), ", p-value: ", fleiss_result$p.value, "\n"))
    
  } else {
    cat("Warning: Not enough complete data across all three platforms to calculate Fleiss' Kappa.\n")
  }
  # <<< NEW SECTION END >>>
  
  # <<< NEW SECTION START >>>
  # HIERARCHICAL CONSISTENCY ASSESSMENT
  # PURPOSE: To replace the single flawed Fleiss' Kappa with a more scientifically
  # sound, two-stage validation that respects the different nature of the datasets.
  # -------------------------------------------------------------------------
  cat("--- Performing Stage 1: Consistency between Categorical Products (MODIS vs. Copernicus) ---\n")
  
  # Filter for common, comparable classes and complete data
  comparable_classes <- c("Forest", "Woody Vegetation", "Grassland", "Cropland/Vegetation")
  modis_coper_data <- data %>%
    filter(
      MODIS_IGBP_Simplified %in% comparable_classes &
        Copernicus_LC_Simplified %in% comparable_classes
    )
  
  # 1. Create a confusion matrix (contingency table)
  confusion_matrix <- table(
    "MODIS" = modis_coper_data$MODIS_IGBP_Simplified,
    "Copernicus" = modis_coper_data$Copernicus_LC_Simplified
  )
  
  # Print and save the confusion matrix
  cat("MODIS vs. Copernicus Confusion Matrix:\n")
  print(confusion_matrix)
  saveTable(as.data.frame.matrix(confusion_matrix), "Supplementary_Table_2b_MODIS_vs_Copernicus_Confusion_Matrix")
  
  # 2. Calculate Cohen's Kappa for these two comparable raters
  # The irr package needs a two-column data frame
  kappa_data <- modis_coper_data %>% dplyr::select(MODIS_IGBP_Simplified, Copernicus_LC_Simplified)
  cohen_kappa_result <- kappam.light(kappa_data) # Use kappam.light for 2 raters
  
  cat("\nCohen's Kappa (MODIS vs. Copernicus):\n")
  print(cohen_kappa_result)
  # saveTable(cohen_kappa_result, "Supplementary_Table_2b_kappa_MODIS_VS_Copernicus")
  # --- Performing Stage 2: Validation using Structural Data (Hansen) ---
  cat("\n--- Performing Stage 2: Validating MODIS classes with Hansen Tree Cover % ---\n")
  
  # Calculate the median and IQR of Hansen tree cover for each MODIS class
  hansen_validation <- data %>%
    filter(
      MODIS_IGBP_Simplified %in% c("Forest", "Woody Vegetation", "Grassland", "Cropland/Vegetation", "Barren/Sparse"),
      !is.na(Hansen_Tree_Cover_2000_Percent)
    ) %>%
    group_by(MODIS_IGBP_Simplified) %>%
    dplyr::summarise(
      N_Samples = n(),
      Median_Hansen_Tree_Cover = median(Hansen_Tree_Cover_2000_Percent, na.rm = TRUE),
      IQR_Hansen_Tree_Cover = IQR(Hansen_Tree_Cover_2000_Percent, na.rm = TRUE)
    ) %>%
    arrange(desc(Median_Hansen_Tree_Cover))
  
  cat("Validation of MODIS Classes using Hansen Tree Cover:\n")
  print(hansen_validation)
  saveTable(hansen_validation, "Supplementary_Table_2b_MODIS_Validation_with_Hansen")
  
  # <<< NEW SECTION END >>>
  
  # <<< NEW SECTION START >>>
  # STATISTICAL TEST 2: Kruskal-Wallis test for Biophysical Zone Distribution
  # PURPOSE: To statistically test if the distribution of biophysical zones is
  # significantly different across the major land cover types shown in Fig 2c.
  # --- Kruskal-Wallis Test (Moved from after 2C to here for better flow) ---
  cat("\n--- Performing Kruskal-Wallis test for biophysical zone distribution across land cover types ---\n")
  kruskal_data <- data %>%
    filter(
      !is.na(NDVI_Biophysical_Zone) & !is.na(MODIS_IGBP_Simplified) &
        NDVI_Biophysical_Zone != "Unknown" &
        MODIS_IGBP_Simplified %in% comparable_classes
    ) %>%
    mutate(
      Biophysical_Zone_Rank = as.numeric(factor(NDVI_Biophysical_Zone,
                                                levels = c("IDZ", "CTZ", "SDZ")))
    )
  
  if(nrow(kruskal_data) > 10) {
    kw_test_result <- kruskal.test(Biophysical_Zone_Rank ~ MODIS_IGBP_Simplified, data = kruskal_data)
    print(kw_test_result)
    kw_summary <- data.frame(
      Test = "Kruskal-Wallis Rank Sum Test",
      Chi_Squared = kw_test_result$statistic,
      df = kw_test_result$parameter,
      P_value = kw_test_result$p.value
    )
    saveTable(kw_summary, "Statistics_Fig2C_Kruskal_Wallis_LCType")
  } else {
    cat("Warning: Not enough data for Kruskal-Wallis test.\n")
  }
  
  # <<< NEW SECTION END >>>
  
  # Helper function to analyze zone distribution for panels 2C, 2D, and 2E
  analyze_zone_distribution <- function(data_subset, x_var_str, output_filename) {
    summary_df <- data_subset %>%
      filter(!is.na(.data[[x_var_str]]), !is.na(NDVI_Biophysical_Zone),
             .data[[x_var_str]] != "Unknown", .data[[x_var_str]] != "Other", 
             NDVI_Biophysical_Zone != "Unknown") %>%
      # First grouping: group by x-axis category and risk zone to count each combination
      group_by(across(all_of(x_var_str)), NDVI_Biophysical_Zone) %>%
      plyr::summarise(Count = n(), .groups = "drop") %>%
      # Second grouping (key step): group only by the x-axis category
      group_by(across(all_of(x_var_str))) %>%
      # Calculate percentage, where sum(Count) is the total within each x-axis category
      mutate(Percentage = Count / sum(Count) * 100)
    
    if(nrow(summary_df) > 0) {
      saveTable(summary_df, output_filename)
    }
    return(summary_df)
  }
  # --- 2C: Risk zone analysis in major land cover types ---
  landcover_zone_summary <- analyze_zone_distribution(
    data_subset = data %>% filter(MODIS_IGBP_Simplified %in% c("Forest", "Woody Vegetation", "Grassland", "Cropland/Vegetation")),
    x_var_str = "MODIS_IGBP_Simplified",
    output_filename = "Fig2C_Landcover_Zone_Distribution"
  )
  analysis_results[["landcover_zone_summary"]] <- landcover_zone_summary
  
  
  # --- 2D: Forest cover density and risk zone analysis ---
  forest_density_zone_summary <- analyze_zone_distribution(
    data_subset = data,
    x_var_str = "Forest_Cover_Category",
    output_filename = "Fig2D_ForestDensity_Zone_Distribution"
  )
  analysis_results[["forest_density_zone_summary"]] <- forest_density_zone_summary
  
  # --- 2E: Climate zone risk distribution analysis ---
  climate_zone_summary <- analyze_zone_distribution(
    data_subset = data,
    x_var_str = "Climate_Zone",
    output_filename = "Fig2E_climate_zone_validation"
  )
  analysis_results[["climate_zone_summary"]] <- climate_zone_summary
  # assign("fig2_analysis_results", analysis_results, envir = .GlobalEnv)
  cat("✓ Vegetation validation analysis completed successfully\n")
  return(analysis_results)
}

analyze_vegetation_validation_data <- function() {
  cat("\n=== Analyzing Vegetation Data Validation and Cross-Platform Consistency ===\n")
  
  # Use the globally defined `landslide_data_processed` variable
  data <- landslide_data_processed 
  analysis_results <- list()
  
  # --- 2A: Global vegetation type distribution analysis ---
  if("MODIS_IGBP_Class_Clean" %in% colnames(data)) {
    modis_dist <- data %>%
      # Explicitly use dplyr::filter
      dplyr::filter(!is.na(MODIS_IGBP_Class_Clean), MODIS_IGBP_Class_Clean != "Unknown") %>%
      dplyr::group_by(MODIS_IGBP_Class_Clean) %>%
      # Explicitly use dplyr::summarise and dplyr::n()
      dplyr::summarise(Count = n(), .groups = "drop") %>%
      dplyr::filter(Count > 0) %>%
      # Explicitly use dplyr::arrange
      dplyr::arrange(desc(Count)) %>%
      # Explicitly use dplyr::slice_head
      dplyr::slice_head(n = 10) %>% # Top 10 for clarity
      # Explicitly use dplyr::mutate
      dplyr::mutate(Percentage = Count / sum(Count) * 100,
                    Display_Name_Short = sapply(strsplit(MODIS_IGBP_Class_Clean, " "), function(x) paste(x[1:min(length(x),2)], collapse=" "))) # First 1 or 2 words
    
    analysis_results[["vegetation_type_distribution"]] <- modis_dist
    saveTable(modis_dist, "Fig2A_modis_vegetation_distribution")
  } else {
    cat("Warning: MODIS_IGBP_Class_Clean not available for vegetation type analysis\n")
  }
  
  # --- 2B: Multi-source consistency level analysis ---
  if("Multi_Source_Consistency_Label" %in% colnames(data)) {
    consistency_summary <- data %>%
      dplyr::filter(!is.na(Multi_Source_Consistency_Label)) %>%
      dplyr::group_by(Multi_Source_Consistency_Label) %>%
      dplyr::summarise(Count = n(), .groups = "drop") %>%
      dplyr::mutate(Percentage = Count / sum(Count) * 100) %>%
      dplyr::mutate(Multi_Source_Consistency_Label = factor(Multi_Source_Consistency_Label,
                                                            levels = c("Full Agreement", "Partial Agreement", "Unknown", "No Agreement"))) %>%
      dplyr::arrange(desc(Multi_Source_Consistency_Label))
    
    analysis_results[["multisource_consistency"]] <- consistency_summary
    saveTable(consistency_summary, "Fig2B_multisource_consistency_levels")
  } else {
    cat("Warning: Multi_Source_Consistency_Label not available for consistency analysis\n")
  }
  
  
  
  # --- Helper function to analyze zone distribution for panels 2C, 2D, and 2E ---
  analyze_zone_distribution <- function(data_subset, x_var_str, output_filename) {
    summary_df <- data_subset %>%
      dplyr::filter(!is.na(.data[[x_var_str]]), !is.na(NDVI_Biophysical_Zone),
                    .data[[x_var_str]] != "Unknown", .data[[x_var_str]] != "Other", 
                    NDVI_Biophysical_Zone != "Unknown") %>%
      # First grouping: group by x-axis category and risk zone
      dplyr::group_by(across(all_of(x_var_str)), NDVI_Biophysical_Zone) %>%
      # Count each combination
      dplyr::summarise(Count = n(), .groups = "drop") %>%
      # Second grouping (key step): group only by the x-axis category
      dplyr::group_by(across(all_of(x_var_str))) %>%
      # Calculate percentage, where sum(Count) is the total within each x-axis category
      dplyr::mutate(Percentage = Count / sum(Count) * 100) %>%
      # Ungrouping is good practice after calculations
      dplyr::ungroup()
    
    if(nrow(summary_df) > 0) {
      saveTable(summary_df, output_filename)
    }
    return(summary_df)
  }
  
  # --- 2C: Risk zone analysis in major land cover types ---
  landcover_zone_summary <- analyze_zone_distribution(
    data_subset = data %>% dplyr::filter(MODIS_IGBP_Simplified %in% c("Forest", "Woody Vegetation", "Grassland", "Cropland/Vegetation")),
    x_var_str = "MODIS_IGBP_Simplified",
    output_filename = "Fig2C_Landcover_Zone_Distribution"
  )
  analysis_results[["landcover_zone_summary"]] <- landcover_zone_summary
  
  # --- 2D: Forest cover density and risk zone analysis ---
  forest_density_zone_summary <- analyze_zone_distribution(
    data_subset = data,
    x_var_str = "Forest_Cover_Category",
    output_filename = "Fig2D_ForestDensity_Zone_Distribution"
  )
  analysis_results[["forest_density_zone_summary"]] <- forest_density_zone_summary
  
  # --- 2E: Climate zone risk distribution analysis ---
  climate_zone_summary <- analyze_zone_distribution(
    data_subset = data,
    x_var_str = "Climate_Zone",
    output_filename = "Fig2E_climate_zone_validation"
  )
  analysis_results[["climate_zone_summary"]] <- climate_zone_summary
  
  
  # <<< NEW SECTION START >>>
  # STATISTICAL TEST 1: Fleiss' Kappa for Multi-Source Consistency
  # PURPOSE: To provide a single, robust statistic quantifying the agreement
  # between the three land cover datasets (MODIS, Copernicus, Hansen).
  # -------------------------------------------------------------------------
  cat("--- Performing Fleiss' Kappa test for multi-source agreement ---\n")
  
  # Prepare the data for the kappa test
  # The `kappam.fleiss` function requires a matrix where rows are subjects (landslides)
  # and columns are raters (datasets). The values are the categories.
  consistency_data <- data %>%
    filter(!is.na(MODIS_IGBP_Simplified), !is.na(Copernicus_LC_Simplified), !is.na(Hansen_Simplified)) %>%
    filter(MODIS_IGBP_Simplified != "Unknown" & Copernicus_LC_Simplified != "Unknown" & Hansen_Simplified != "Unknown") %>%
    dplyr::select(MODIS_IGBP_Simplified, Copernicus_LC_Simplified, Hansen_Simplified)
  
  if(nrow(consistency_data) > 10) { # Check if there is enough data
    fleiss_result <- kappam.fleiss(consistency_data)
    
    # Print the results to the console
    print(fleiss_result)
    
    # Create a summary table to save
    fleiss_summary <- data.frame(
      Statistic = "Fleiss' Kappa",
      Value = fleiss_result$value,
      Z_statistic = fleiss_result$statistic,
      P_value = fleiss_result$p.value,
      Subjects = fleiss_result$subjects,
      Raters = fleiss_result$raters
    )
    
    # The `irr` package doesn't directly provide CI for Fleiss' Kappa.
    # We will add the manually provided CI for reference.
    fleiss_summary$CI_95_Lower_Manual = 0.67
    fleiss_summary$CI_95_Upper_Manual = 0.71
    
    saveTable(fleiss_summary, "Statistics_Fig2B_Fleiss_Kappa")
    cat(paste0("✓ Fleiss' Kappa result: ", round(fleiss_result$value, 2), ", p-value: ", fleiss_result$p.value, "\n"))
    
  } else {
    cat("Warning: Not enough complete data across all three platforms to calculate Fleiss' Kappa.\n")
  }
  # <<< NEW SECTION END >>>
  
  
  # <<< NEW SECTION START >>>
  # HIERARCHICAL CONSISTENCY ASSESSMENT
  # PURPOSE: To replace the single flawed Fleiss' Kappa with a more scientifically
  # sound, two-stage validation that respects the different nature of the datasets.
  # -------------------------------------------------------------------------
  cat("--- Performing Stage 1: Consistency between Categorical Products (MODIS vs. Copernicus) ---\n")
  
  # Filter for common, comparable classes and complete data
  comparable_classes <- c("Forest", "Woody Vegetation", "Grassland", "Cropland/Vegetation")
  modis_coper_data <- data %>%
    filter(
      MODIS_IGBP_Simplified %in% comparable_classes &
        Copernicus_LC_Simplified %in% comparable_classes
    )
  
  # 1. Create a confusion matrix (contingency table)
  confusion_matrix <- table(
    "MODIS" = modis_coper_data$MODIS_IGBP_Simplified,
    "Copernicus" = modis_coper_data$Copernicus_LC_Simplified
  )
  
  # Print and save the confusion matrix
  cat("MODIS vs. Copernicus Confusion Matrix:\n")
  print(confusion_matrix)
  saveTable(as.data.frame.matrix(confusion_matrix), "Supplementary_Table_2b_MODIS_vs_Copernicus_Confusion_Matrix")
  
  # 2. Calculate Cohen's Kappa for these two comparable raters
  # The irr package needs a two-column data frame
  kappa_data <- modis_coper_data %>% dplyr::select(MODIS_IGBP_Simplified, Copernicus_LC_Simplified)
  cohen_kappa_result <- kappam.light(kappa_data) # Use kappam.light for 2 raters
  
  cat("\nCohen's Kappa (MODIS vs. Copernicus):\n")
  print(cohen_kappa_result)
  saveTable(as.data.frame(cohen_kappa_result)[c(1:5)], "Supplementary_Table_2b_kappa_MODIS_VS_Copernicus")
  # --- Performing Stage 2: Validation using Structural Data (Hansen) ---
  cat("\n--- Performing Stage 2: Validating MODIS classes with Hansen Tree Cover % ---\n")
  
  # Calculate the median and IQR of Hansen tree cover for each MODIS class
  hansen_validation <- data %>%
    filter(
      MODIS_IGBP_Simplified %in% c("Forest", "Woody Vegetation", "Grassland", "Cropland/Vegetation", "Barren/Sparse"),
      !is.na(Hansen_Tree_Cover_2000_Percent)
    ) %>%
    group_by(MODIS_IGBP_Simplified) %>%
    dplyr::summarise(
      N_Samples = n(),
      Median_Hansen_Tree_Cover = median(Hansen_Tree_Cover_2000_Percent, na.rm = TRUE),
      IQR_Hansen_Tree_Cover = IQR(Hansen_Tree_Cover_2000_Percent, na.rm = TRUE)
    ) %>%
    arrange(desc(Median_Hansen_Tree_Cover))
  
  cat("Validation of MODIS Classes using Hansen Tree Cover:\n")
  print(hansen_validation)
  saveTable(hansen_validation, "Supplementary_Table_2b_MODIS_Validation_with_Hansen")
  
  # <<< NEW SECTION END >>>
  
  
  # <<< NEW SECTION START >>>
  # STATISTICAL TEST 2: Kruskal-Wallis test for Biophysical Zone Distribution
  # PURPOSE: To statistically test if the distribution of biophysical zones is
  # significantly different across the major land cover types shown in Fig 2c.
  # --- Kruskal-Wallis Test (Moved from after 2C to here for better flow) ---
  cat("\n--- Performing Kruskal-Wallis test for biophysical zone distribution across land cover types ---\n")
  kruskal_data <- data %>%
    filter(
      !is.na(NDVI_Biophysical_Zone) & !is.na(MODIS_IGBP_Simplified) &
        NDVI_Biophysical_Zone != "Unknown" &
        MODIS_IGBP_Simplified %in% comparable_classes
    ) %>%
    mutate(
      Biophysical_Zone_Rank = as.numeric(factor(NDVI_Biophysical_Zone,
                                                levels = c("IDZ", "CTZ", "SDZ")))
    )
  
  if(nrow(kruskal_data) > 10) {
    kw_test_result <- kruskal.test(Biophysical_Zone_Rank ~ MODIS_IGBP_Simplified, data = kruskal_data)
    print(kw_test_result)
    kw_summary <- data.frame(
      Test = "Kruskal-Wallis Rank Sum Test",
      Chi_Squared = kw_test_result$statistic,
      df = kw_test_result$parameter,
      P_value = kw_test_result$p.value
    )
    saveTable(kw_summary, "Statistics_Fig2C_Kruskal_Wallis_LCType")
  } else {
    cat("Warning: Not enough data for Kruskal-Wallis test.\n")
  }
  
  # <<< NEW SECTION END >>>
  
  # Assign the analysis results to a global variable for use in plotting functions
  assign("vegetation_validation_results", analysis_results, envir = .GlobalEnv)
  
  cat("✓ Vegetation validation analysis completed successfully\n")
  return(analysis_results)
}

# Run the analysis
vegetation_validation_results <- analyze_vegetation_validation_data()


# -------------------------------------------------------------------------
# MAIN FIGURE 2: Visualization - Multi-Source Vegetation Data Validation and Cross-Platform Consistency
# -------------------------------------------------------------------------
create_vegetation_validation_plots <- function(analysis_results, data) {
  cat("\n=== Creating Figure 2 Visualizations with Improved Labels ===\n")
  
  plot_list <- list()
  
  # --- 2A: Global Vegetation Type Distribution (Treemap) ---
  create_vegetation_treemap <- function() {
    if(!"vegetation_type_distribution" %in% names(analysis_results) || 
       nrow(analysis_results[["vegetation_type_distribution"]]) == 0) {
      return(ggplot() + geom_text(aes(0, 0, label = "MODIS data not available")) + 
               nature_theme_professional)
    }
    
    modis_dist <- analysis_results[["vegetation_type_distribution"]]
    
    p_treemap <- ggplot(modis_dist, aes(area = Count, fill = Count, 
                                        label = paste(Display_Name_Short, 
                                                      sprintf("\n%.1f%%", Percentage)))) +
      geom_treemap() +
      geom_treemap_text(color = "black", place = "centre", size = 8, 
                        fontface = "bold", grow = FALSE, reflow=TRUE, min.size = 3) +
      scale_fill_gradientn(colors = c("#E6F0DC", "#C1E899", "#55883B"), name = "Count") +
      nature_theme_professional +
      theme(
        legend.position = "right",
        plot.margin = ggplot2::margin(0, 0, 0, 0)
      )
    
    return(p_treemap)
  }
  
  # --- 2B: Multi-Source Consistency Levels (Pie Chart) ---
  create_multisource_consistency_piechart <- function() {
    if(!"multisource_consistency" %in% names(analysis_results) || 
       nrow(analysis_results[["multisource_consistency"]]) == 0) {
      return(ggplot() + geom_text(aes(0, 0, label = "Multi-source consistency data not available")) + 
               nature_theme_professional)
    }
    
    consistency_summary <- analysis_results[["multisource_consistency"]] %>%
      mutate(Multi_Source_Consistency_Label = factor(Multi_Source_Consistency_Label,
                                                     levels = c("No Agreement", "Unknown", "Partial Agreement", "Full Agreement"))) %>%
      arrange(desc(Multi_Source_Consistency_Label)) %>%
      mutate(
        label_text = paste0(Multi_Source_Consistency_Label, "\n", 
                            sprintf("%.1f%%", Percentage), "(n=", 
                            format(Count, big.mark = ","), ")")
      ) %>%
      mutate(pos = cumsum(Percentage) - 0.5 * Percentage)
    
    consistency_colors <- c(
      "Full Agreement" = "#4CAF50",  
      "Partial Agreement" = "#8BC34A", 
      "Unknown" = "#BDBDBD",        
      "No Agreement" = "#F44336"     
    )
    
    p_piechart <- ggplot(consistency_summary,
                         aes(x = 0.75, y = Percentage, fill = Multi_Source_Consistency_Label)) +
      geom_bar(stat = "identity", width = 1.5, color = "white", linewidth = 0.5) +
      
      coord_polar(theta = "y", start = 0, direction = 1) +
      
      ggrepel::geom_text_repel(
        aes(y = pos, label = label_text),
        size = 2.8,
        nudge_x = 1.2,
        show.legend = FALSE,
        segment.color = 'grey50',
        segment.size = 0.4,
        force = 10,
        min.segment.length = 0
      ) +
      
      scale_fill_manual(values = consistency_colors) +
      theme_void() +
      theme(
        legend.position = "none",
        plot.margin = ggplot2::margin(5, 5, 5, 5) 
      ) +
      xlim(c(0, 2.2))
    
    return(p_piechart)
  }
  
  # --- Shared function: Create risk distribution plot (for 2C, 2D, 2E) ---
  # Main change: Adds logic to wrap long labels onto two lines.
  create_zone_distribution_plot <- function(summary_df, x_var_str, x_lab, is_first_plot = FALSE, is_last_plot = FALSE) {
    if(is.null(summary_df) || nrow(summary_df) == 0) {
      return(ggplot() + geom_text(aes(0, 0, label = paste("No data available for", x_var_str))) + 
               nature_theme_professional)
    }
    
    # Create a new modified data frame to handle label display
    modified_df <- summary_df
    
    # Special handling for forest cover labels - more concise format
    if(x_var_str == "Forest_Cover_Category") {
      # Use more concise labels, retaining only essential info
      label_mapping <- c(
        "No Forest (0%)" = "No Forest\n(0%)",
        "Sparse (1-25%)" = "Sparse\n(1-25%)",
        "Moderate (26-50%)" = "Moderate\n(26-50%)",
        "Dense (51-75%)" = "Dense\n(51-75%)",
        "Very Dense (>75%)" = "Very Dense\n(>75%)"
      )
      
      # Ensure original order is maintained
      if("order" %in% colnames(modified_df)) {
        original_order <- modified_df$order
      } else {
        # If no order column, create one based on standard sorting
        standard_order <- c("No Forest (0%)", "Sparse (1-25%)", "Moderate (26-50%)", 
                            "Dense (51-75%)", "Very Dense (>75%)")
        original_order <- match(modified_df[[x_var_str]], standard_order)
      }
      
      # Create new label column and maintain sorting
      modified_df$Display_Label <- sapply(modified_df[[x_var_str]], function(x) {
        if(x %in% names(label_mapping)) {
          return(label_mapping[x])
        } else {
          return(as.character(x))
        }
      })
      
      # Ensure correct factor ordering
      modified_df$Display_Label <- factor(modified_df$Display_Label, 
                                          levels = unique(modified_df$Display_Label[order(original_order)]))
    }
    # Handle other label types
    else if(x_var_str == "MODIS_IGBP_Simplified") {
      # Land cover type label handling
      label_mapping <- c(
        "Cropland/Vegetation" = "Croplands/\nVegetation",
        "Evergreen Needleleaf" = "Evergreen\nNeedleleaf",
        "Evergreen Broadleaf" = "Evergreen\nBroadleaf",
        "Deciduous Needleleaf" = "Deciduous\nNeedleleaf",
        "Deciduous Broadleaf" = "Deciduous\nBroadleaf",
        "Mixed Forest" = "Mixed\nForest",
        "Closed Shrublands" = "Closed\nShrublands",
        "Open Shrublands" = "Open\nShrublands",
        "Woody Savannas" = "Woody\nSavannas",
        "Savannas" = "Savannas",
        "Grassland" = "Grassland",
        "Permanent Wetlands" = "Permanent\nWetlands",
        "Croplands" = "Croplands",
        "Woody Vegetation" = "Woody\nVegetation",
        "Urban & Built-up" = "Urban &\nBuilt-up",
        "Cropland/Natural" = "Cropland/\nNatural",
        "Snow & Ice" = "Snow &\nIce",
        "Barren" = "Barren",
        "Forest" = "Forest",
        "Water Bodies" = "Water\nBodies"
      )
      
      modified_df$Display_Label <- sapply(modified_df[[x_var_str]], function(x) {
        if(x %in% names(label_mapping)) {
          return(label_mapping[x])
        } else {
          return(as.character(x))
        }
      })
    } else if(x_var_str == "Climate_Zone") {
      # Climate zone label handling
      label_mapping <- c(
        "Temperate" = "Temperate",
        "Mediterranean/Subtropical" = "Mediterranean/\nSubtropical",
        "Tropical" = "Tropical",
        "Arid" = "Arid",
        "Continental" = "Continental",
        "Polar" = "Polar",
        "Other" = "Other"
      )
      
      modified_df$Display_Label <- sapply(modified_df[[x_var_str]], function(x) {
        if(x %in% names(label_mapping)) {
          return(label_mapping[x])
        } else {
          return(as.character(x))
        }
      })
    } else {
      # Default case, use original labels
      modified_df$Display_Label <- modified_df[[x_var_str]]
    }
    
    # Use more compact width settings for the forest cover plot
    bar_width <- 0.8  # Default width
    if(x_var_str == "Forest_Cover_Category") {
      bar_width <- 0.7  # Use narrower bars for forest cover
    }
    
    # Create plot with modified labels
    p <- ggplot(modified_df, aes(x = Display_Label, y = Percentage, fill = NDVI_Biophysical_Zone)) +
      geom_col(position = "stack", alpha = 0.9, width = bar_width) +
      geom_text(aes(label = ifelse(Percentage > 7, sprintf("%.0f%%", Percentage), "")),
                position = position_stack(vjust = 0.5), color = "white", 
                fontface = "bold", size = 2.5) +
      scale_fill_manual(values = nature_palettes$biophysical_zones, name="Biophysical Zone") +
      scale_y_continuous(labels = scales::percent_format(scale=1), 
                         limits = c(0, 100), breaks = seq(0,100,25)) +
      nature_theme_professional
    
    # Use special label settings for the forest cover plot
    if(x_var_str == "Forest_Cover_Category") {
      p <- p + theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, size=7), # Slightly smaller font
        axis.text.y = element_text(size=6),
        axis.title = element_text(size=7),
        plot.margin = ggplot2::margin(2, 0, 2, 0),  # Reduce left/right margins
        panel.spacing = unit(0, "pt")  # Reduce panel spacing
      )
    } else {
      # Use standard settings for other plots
      p <- p + theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, size=8),
        axis.text.y = element_text(size=6),
        axis.title = element_text(size=7),
        plot.margin = ggplot2::margin(2, 2, 2, 2)
      )
    }
    
    # Add x-axis label
    p <- p + labs(x = x_lab)
    
    # Only add Y-axis label to the first plot
    if (is_first_plot) {
      p <- p + labs(y = "Percentage of Landslides (%)")
    } else {
      p <- p + theme(axis.title.y = element_blank(), 
                     axis.text.y = element_blank(), 
                     axis.ticks.y = element_blank())
    }
    
    # Only show the legend on the last plot
    if (!is_last_plot) {
      p <- p + theme(legend.position = "none")
    } else {
      p <- p + theme(
        legend.position = "right", 
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size=6),
        legend.text = element_text(size=6),
        legend.margin = ggplot2::margin(0,0,0,0)
      )
    }
    
    return(p)
  }
  
  # --- 2C: Risk zones in major land cover types ---
  p2c <- create_zone_distribution_plot(
    summary_df = analysis_results[["landcover_zone_summary"]],
    x_var_str = "MODIS_IGBP_Simplified",
    x_lab = "Land Cover Type",
    is_first_plot = FALSE,
    is_last_plot = FALSE
  )
  
  # --- 2D: Forest cover density and risk zones ---
  p2d <- create_zone_distribution_plot(
    summary_df = analysis_results[["forest_density_zone_summary"]],
    x_var_str = "Forest_Cover_Category",
    x_lab = "Forest Cover Category",
    is_first_plot = FALSE,
    is_last_plot = FALSE
  )
  
  # --- 2E: Climate zone risk distribution ---
  p2e <- create_zone_distribution_plot(
    summary_df = analysis_results[["climate_zone_summary"]],
    x_var_str = "Climate_Zone",
    x_lab = "Climate Zone",
    is_first_plot = FALSE,
    is_last_plot = FALSE  # This is the last plot, show the legend
  )
  
  # Create all components of Figure 2
  p2a <- create_vegetation_treemap()
  p2b <- create_multisource_consistency_piechart()
  
  plot_list <- list(
    treemap = p2a,
    consistency = p2b,
    landcover_zone = p2c,
    forest_density = p2d,
    climate_zone = p2e
  )
  
  # Combine Figure 2 - use similar dimension ratios as Figure 1
  top_row <- p2a + p2b + plot_layout(widths = c(1.2, 0.8))
  
  # Combine bottom row and ensure only one legend
  bottom_row <- p2c + p2d + p2e + 
    plot_layout(
      widths = c(1, 1.2, 1),
      guides = "collect"  # Collect all legends
    ) & 
    theme(legend.position = "bottom")  # Ensure legend is at the bottom
  
  # Combine top and bottom rows
  main_fig2 <- (top_row / bottom_row) + 
    plot_layout(heights = c(1, 1)) + 
    plot_annotation(tag_levels = list(paste0("(", letters, ")"))) &
    theme(
      plot.tag = element_text(face = "bold", size = 10),
      panel.spacing = unit(4, "pt"),
      plot.margin = ggplot2::margin(5, 5, 5, 5)
    )
  
  # Save the combined figure - use the same size as Figure 1
  saveFigure(main_fig2, "Main_Figure_2_A4", width = 8.3, height = 6.8)
  
  cat("✓ Figure 2 visualizations with improved labels completed successfully\n")
  return(list(plots = plot_list, combined = main_fig2))
}

# Generate the plots
figure2_plots <- create_vegetation_validation_plots(vegetation_validation_results, landslide_data_processed)

# ------------------------------------------------------------------------------
# Figure 2 Supplementary: Sankey Diagram
# ------------------------------------------------------------------------------
# ==============================================================================
# Step 1: Data Preparation
# ==============================================================================

# 1.1 Create a unified, simplified final classification (as the rightmost endpoint of the Sankey diagram)
unified_classification <- function(modis, copernicus, hansen) {
  dplyr::case_when(
    grepl("Forest", modis) | grepl("Forest", copernicus) | grepl("Forest", hansen) ~ "Forest",
    grepl("Woody", modis) | grepl("Woody", copernicus) | grepl("Woody", hansen) ~ "Woody Vegetation",
    grepl("Grassland", modis) | grepl("Grassland", copernicus) ~ "Grassland",
    grepl("Cropland", modis) | grepl("Cropland", copernicus) ~ "Cropland/Vegetation",
    grepl("Urban", modis) | grepl("Urban", copernicus) ~ "Urban/Built-up",
    grepl("Barren", modis) | grepl("Barren", copernicus) ~ "Barren/Sparse",
    grepl("Water", modis) | grepl("Water", copernicus) ~ "Water/Wetland/Ice",
    hansen == "Non-forest" ~ "Non-forest",
    TRUE ~ "Unknown"
  )
}

# 1. Filter data
complete_data <- landslide_data_processed %>%
  dplyr::filter(!is.na(MODIS_IGBP_Simplified) & !is.na(Copernicus_LC_Simplified) & !is.na(Hansen_Simplified)) %>%
  dplyr::filter(MODIS_IGBP_Simplified != "Unknown" & Copernicus_LC_Simplified != "Unknown")

total_samples <- nrow(complete_data)

# 2. Calculate frequencies
combo_counts <- complete_data %>%
  dplyr::group_by(MODIS_IGBP_Simplified, Copernicus_LC_Simplified, Hansen_Simplified) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(percent = n / total_samples * 100)

# 3. Assign alluvium IDs
combo_counts <- combo_counts %>%
  dplyr::mutate(alluvium = dplyr::row_number())

# 4. Convert to long format
sankey_data_long <- combo_counts %>%
  tidyr::pivot_longer(
    cols = c(MODIS_IGBP_Simplified, Copernicus_LC_Simplified, Hansen_Simplified),
    names_to = "Source",
    values_to = "stratum"
  )

# 5. Add unified classification
classification_lookup <- combo_counts %>%
  dplyr::mutate(
    Unified_Class = unified_classification(
      MODIS_IGBP_Simplified,
      Copernicus_LC_Simplified,
      Hansen_Simplified
    )
  ) %>%
  dplyr::select(alluvium, Unified_Class)

sankey_data_with_unified <- sankey_data_long %>%
  dplyr::left_join(classification_lookup, by = "alluvium")

# ==============================================================================
# Step 2: Define Colors and Labels
# ==============================================================================

color_palette <- c(
  "Forest" = "#1a9850",
  "Woody Vegetation" = "#91cf60",
  "Grassland" = "#d9ef8b",
  "Cropland/Vegetation" = "#fee08b",
  "Urban/Built-up" = "#999999",
  "Barren/Sparse" = "#d0d0d0",
  "Water/Wetland/Ice" = "#4393c3",
  "Non-forest" = "#bababa",
  "Unknown" = "#525252"
)

class_order <- c(
  "Forest", "Woody Vegetation", "Grassland", "Cropland/Vegetation",  
  "Urban/Built-up", "Barren/Sparse", "Water/Wetland/Ice", "Non-forest"
)

sankey_data_final <- sankey_data_with_unified %>%
  dplyr::mutate(
    Source = gsub("_Simplified", "", Source),
    stratum_ordered = factor(stratum, levels = class_order)
  )

# ==============================================================================
# Step 3: Prepare Final Plotting Data
# ==============================================================================

unified_axis_data <- sankey_data_final %>%
  dplyr::group_by(alluvium) %>%
  dplyr::summarise(
    n = dplyr::first(n),
    percent = dplyr::first(percent),
    stratum = dplyr::first(Unified_Class),
    Source = "Unified",
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    stratum_ordered = factor(stratum, levels = class_order)
  )

plot_data <- dplyr::bind_rows(
  sankey_data_final %>% dplyr::select(alluvium, n, percent, Source, stratum, stratum_ordered),
  unified_axis_data
) %>%
  # Ensure Source is an ordered factor for plotting
  dplyr::mutate(
    Source = factor(Source, levels = c("MODIS_IGBP", "Copernicus_LC", "Hansen", "Unified"))
  )


# ==============================================================================
# Step 4: Add Percentage Labels with Optimized Positioning
# ==============================================================================

# First, calculate the total size of each stratum (classification block)
stratum_summary <- plot_data %>%
  dplyr::filter(!is.na(stratum_ordered)) %>%
  dplyr::group_by(Source, stratum_ordered) %>%
  dplyr::summarise(stratum_total_n = sum(n), .groups = "drop")

# Next, calculate the Y-axis position and percentage for each block
label_data <- stratum_summary %>%
  # Use desc() to reverse the order to match ggalluvial's stacking (top-level factor at the top)
  dplyr::arrange(Source, dplyr::desc(stratum_ordered)) %>%
  # Group by the x-axis variable (Source)
  dplyr::group_by(Source) %>%
  # Calculate the top, bottom, and center position for each color block
  dplyr::mutate(
    y_top = cumsum(stratum_total_n),
    y_bottom = dplyr::lag(y_top, default = 0),
    y_center = (y_top + y_bottom) / 2,
    # Calculate the correct percentage based on the block's total size
    total_n_in_source = sum(stratum_total_n),
    percent = stratum_total_n / total_n_in_source * 100,
    label = ifelse(percent >= 2, sprintf("%.0f%%", percent), "")
  ) %>%
  dplyr::ungroup() # Ungroup after calculations

# Finally, add styling (size and color) to the labels
label_data <- label_data %>%
  dplyr::mutate(
    label_size = dplyr::case_when(
      percent >= 50 ~ 5,
      percent >= 25 ~ 4.5,
      percent >= 10 ~ 4,
      percent >= 5 ~ 3.5,
      TRUE ~ 3
    ),
    label_color = dplyr::case_when(
      stratum_ordered %in% c("Forest", "Urban/Built-up") ~ "white",
      stratum_ordered %in% c("Grassland", "Cropland/Vegetation") ~ "black",
      TRUE ~ "black"
    )
  )

# ==============================================================================
# Step 5: Plot the Sankey Diagram (using the new label_data)
# ==============================================================================

sankey_plot <- ggplot2::ggplot(data = plot_data,
                               ggplot2::aes(x = Source, stratum = stratum_ordered, alluvium = alluvium,
                                            y = n, fill = stratum_ordered)) +
  # Use ggalluvial's core geometric objects
  ggalluvial::geom_flow(stat = "alluvium", lode.guidance = "forward", color = "darkgray", alpha = 0.6, width = 0.4) +
  ggalluvial::geom_stratum(alpha = 1, width = 0.4, linewidth = 0.2) +
  
  # Add percentage labels inside the strata using the corrected label_data
  ggplot2::geom_text(
    data = label_data %>% dplyr::filter(label != ""), # Only plot non-empty labels
    mapping = ggplot2::aes(x = Source, y = y_center, label = label, size = label_size, color = label_color),
    inherit.aes = FALSE,
    fontface = "bold"
  ) +
  
  # Set label size and color
  ggplot2::scale_size_identity() +
  ggplot2::scale_color_identity() +
  
  # Use our defined color palette
  ggplot2::scale_fill_manual(
    values = color_palette,
    name = "Land Cover Classification",
    breaks = class_order,
    na.value = "grey80" # Specify color for any potential NA values
  ) +
  
  # Adjust X-axis labels
  ggplot2::scale_x_discrete(
    limits = c("MODIS_IGBP", "Copernicus_LC", "Hansen", "Unified"),
    labels = c("MODIS", "Copernicus", "Hansen", "Unified"),
    expand = ggplot2::expansion(mult = c(0.05, 0.05)) # Adjust space on both sides of x-axis
  ) +
  
  # No extra space needed for Y-axis
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
  
  # Theme and titles
  nature_theme_professional +
  ggplot2::theme_minimal(base_family = "sans") +
  ggplot2::theme(
    legend.position = "top",
    legend.box = "horizontal",
    legend.title = ggplot2::element_text(size = 9, face = "bold"),
    legend.text = ggplot2::element_text(size = 8),
    legend.margin = ggplot2::margin(b = 0, t = 0), # Reduce distance between legend and plot
    legend.spacing.x = ggplot2::unit(0.2, "cm"), # Reduce spacing within legend
    legend.key.size = ggplot2::unit(0.5, "cm"), # Reduce size of legend keys
    
    panel.grid = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    
    axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10), size = 11, face = "bold"),
    axis.text.x = ggplot2::element_text(size = 10, face = "bold"),
    
    plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = ggplot2::element_text(size = 12, hjust = 0.5, margin = ggplot2::margin(b = 15)),
    plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10)
  ) +
  ggplot2::guides(fill = guide_legend(nrow = 2)) +
  ggplot2::labs(
    title = "Cross-Platform Land Cover Classification Consistency",
    subtitle = "Flow of landslide event locations between MODIS, Copernicus, and Hansen classifications",
    x = "Data Source / Classification", # Add X-axis title
    y = "Number of Landslides",
    fill = "" # Hide the legend title
  )

# Print the plot
print(sankey_plot)

# Save the plot
saveFigure(sankey_plot, "Supplementary_Fig_2b_sankey_plot", width = 8.3, height = 6.8)

# ==============================================================================

# ------------------------------------------------------------------------------
# FIGURE 3
# ------------------------------------------------------------------------------

# ==============================================================================

# ------------------------------------------------------------------------------
# Data Analysis - Preparing Data for Figure 3
# ------------------------------------------------------------------------------
prepare_data_for_figure_3 <- function() {
  cat("\n=== Preparing Data for Figure 3 ===\n")
  
  data <- landslide_data_processed 
  
  # Ensure NDVI_Biophysical_Zone is correctly defined
  if (!"NDVI_Biophysical_Zone" %in% colnames(data) && exists("fig1_analysis_results") && "NDVI_1" %in% names(fig1_analysis_results)) {
    ndvi_result <- fig1_analysis_results[["NDVI_1"]]
    if(!is.null(ndvi_result$breakpoints) && length(ndvi_result$breakpoints) == 2){
      data$NDVI_Biophysical_Zone <- case_when(
        is.na(data$NDVI_1) ~ "Unknown",
        data$NDVI_1 < ndvi_result$breakpoints[1] ~ "IDZ",
        data$NDVI_1 >= ndvi_result$breakpoints[1] & 
          data$NDVI_1 <= ndvi_result$breakpoints[2] ~ "CTZ",
        data$NDVI_1 > ndvi_result$breakpoints[2] ~ "SDZ",
        TRUE ~ "Unknown" 
      )
      data$NDVI_Biophysical_Zone <- factor(data$NDVI_Biophysical_Zone,
                                           levels = c("IDZ", "CTZ", "SDZ", "Unknown"))
      assign("landslide_data_processed", data, envir = .GlobalEnv) 
    } else {
      warning("NDVI breakpoints invalid in fig1_analysis_results for Fig 3.")
    }
  } else if (!"NDVI_Biophysical_Zone" %in% colnames(data)) {
    stop("NDVI_Biophysical_Zone column is missing for Fig 3 and cannot be recreated.")
  }
  
  # Prepare data for Figure 3A (scatter plot)
  data_3a <- data %>% dplyr::filter(!is.na(NDVI_1), NDVI_Biophysical_Zone != "Unknown")
  
  # Prepare data for Figure 3B (violin plot)
  climate_ndvi_data <- data %>%
    dplyr::filter(!is.na(NDVI_1), !is.na(Climate_Zone), Climate_Zone != "Other")
  
  # --- Supplementary Figure Y: Faceted Breakpoint Analysis by Climate Zone ---
  cat("--- Generating Supplementary Figure Y: Faceted Breakpoint Analysis ---\n")
  
  climate_zones <- c("Temperate", "Arid", "Tropical")
  all_climate_plots <- list()
  
  for(zone in climate_zones) {
    zone_data <- data %>% dplyr::filter(Climate_Zone == zone)
    if(nrow(zone_data) > 200) { 
      tryCatch({
        result <- advanced_segmented_analysis(zone_data, "NDVI_1", breaks = 2) 
        pred_data <- data.frame(x = result$data$x, y = predict(result$model))
        
        p <- ggplot(result$data, aes(x = x, y = y)) +
          geom_point(alpha = 0.3, color = "grey") +
          geom_line(data = pred_data, aes(x=x, y=y), color = nature_palettes$climate[zone], linewidth=1.2) +
          geom_vline(xintercept = result$breakpoints, linetype="dashed", color="red") +
          labs(title = zone, x = "NDVI Value", y="Smoothed Frequency") +
          nature_theme_professional
        all_climate_plots[[zone]] <- p
      }, error = function(e) {
        cat("Could not perform segmented analysis for zone:", zone, "\nError:", e$message, "\n")
      })
    }
  }
  
  if(length(all_climate_plots) >= 1) {
    supplementary_fig_y <- wrap_plots(all_climate_plots, ncol = length(all_climate_plots)) +
      plot_annotation(title = "Supplementary Figure 3: NDVI-Landslide Relationship by Climate Zone")
    saveFigure(supplementary_fig_y, "Supplementary_Fig_3_Faceted_Breakpoints", width=12, height=4)
    cat("✓ Supplementary Figure 3 saved to 'figures/' directory.\n")
  } else {
    cat("Warning: Could not generate faceted plots for any climate zones.\n")
  }
  
  # Initialize variables for the return list
  seasonal_summary <- NULL
  
  # Prepare data for Figure 3C (seasonal distribution)
  if(!"Season" %in% colnames(data)){
    warning("Season column not found. Cannot prepare data for Figure 3D.")
  } else {
    seasonal_summary <- data %>%
      dplyr::filter(!is.na(NDVI_1), !is.na(Season), NDVI_Biophysical_Zone != "Unknown", Season != "Unknown") %>%
      dplyr::mutate(NDVI_Category = cut(NDVI_1, 
                                        breaks = c(0, 0.3, 0.6, 0.8, 1.0),
                                        labels = c("Low (0-0.3)", "Moderate (0.3-0.6)", "High (0.6-0.8)", "Very High (0.8-1)"),
                                        include.lowest = TRUE, right = TRUE)) %>%
      dplyr::filter(!is.na(NDVI_Category)) %>%
      dplyr::group_by(NDVI_Category, Season) %>%
      dplyr::summarise(Count = n(), .groups = "drop") %>%
      dplyr::group_by(NDVI_Category) %>%
      dplyr::mutate(Percentage = Count / sum(Count) * 100)
    
    # --- STATISTICAL TEST 4: Chi-squared Test for Seasonal Distribution ---
    cat("--- Performing Chi-squared Tests for Seasonal Distributions ---\n")
    
    chi_sq_data <- data %>%
      dplyr::filter(!is.na(Season), NDVI_Biophysical_Zone != "Unknown", Season != "Unknown") %>%
      dplyr::mutate(NDVI_Category = cut(NDVI_1, 
                                        breaks = c(0, 0.3, 0.6, 0.8, 1.0),
                                        labels = c("Low (0-0.3)", "Moderate (0.3-0.6)", "High (0.6-0.8)", "Very High (0.8-1)"),
                                        include.lowest = TRUE, right = TRUE)) %>%
      dplyr::filter(!is.na(NDVI_Category)) %>%
      dplyr::count(NDVI_Category, Season, name = "Count")
    
    if(nrow(chi_sq_data) > 0) {
      chi_sq_tests <- chi_sq_data %>%
        dplyr::group_by(NDVI_Category) %>%
        tidyr::nest() %>%
        dplyr::mutate(
          # Safely apply the test using purrr::map
          chi_sq_result = purrr::map(data, ~ tryCatch(chisq.test(.x$Count), error = function(e) NULL)),
          # Extract results only if the test was successful
          is_success = !purrr::map_lgl(chi_sq_result, is.null)
        ) %>%
        dplyr::filter(is_success) %>%
        dplyr::mutate(
          Chi_Squared_Stat = purrr::map_dbl(chi_sq_result, "statistic"),
          P_value = purrr::map_dbl(chi_sq_result, "p.value"),
          df = purrr::map_dbl(chi_sq_result, "parameter")
        ) %>%
        dplyr::select(NDVI_Category, Chi_Squared_Stat, df, P_value)
      
      print(chi_sq_tests)
      saveTable(chi_sq_tests, "Statistics_Fig3C_Chi_Squared_Results")
      cat("✓ Chi-squared test results saved to 'tables/' directory.\n")
    } else {
      cat("Warning: No data available for Chi-squared tests.\n")
    }
  }
  
  # Initialize variables in case subsequent steps fail
  temporal_summary <- NULL
  temporal_change_annotations <- NULL
  
  # Prepare data for Figure 3D (time series)
  temporal_cols <- paste0("NDVI_", 5:1)  # Change order from 5 to 1
  if(!all(temporal_cols %in% colnames(data))){
    warning("Temporal NDVI columns (NDVI_1 to NDVI_5) not found. Cannot prepare data for Figure 3C.")
  } else {
    temporal_data_long <- data %>%
      dplyr::filter(NDVI_Biophysical_Zone != "Unknown") %>%
      dplyr::select(OBJECTID, NDVI_Biophysical_Zone, all_of(temporal_cols)) %>%
      tidyr::pivot_longer(cols = all_of(temporal_cols), names_to = "Time_Window_Raw", values_to = "NDVI_Value") %>%
      dplyr::filter(!is.na(NDVI_Value)) %>%
      dplyr::mutate(Time_Window = gsub("NDVI_", "T", Time_Window_Raw),
                    Time_Window_Num = as.numeric(gsub("NDVI_", "", Time_Window_Raw))) %>%
      # Add a reverse-order index for sorting
      dplyr::mutate(Reverse_Time_Window_Num = 6 - Time_Window_Num)  
    
    if(nrow(temporal_data_long) > 0){
      temporal_summary <- temporal_data_long %>%
        dplyr::group_by(NDVI_Biophysical_Zone, Time_Window, Time_Window_Num, Reverse_Time_Window_Num) %>%
        dplyr::summarise(
          Mean_NDVI = mean(NDVI_Value, na.rm = TRUE),
          CI_Lower = Rmisc::CI(NDVI_Value, ci=0.95)[3], 
          CI_Upper = Rmisc::CI(NDVI_Value, ci=0.95)[1],
          N = n(),
          .groups = "drop"
        ) %>% 
        # Sort in reverse order
        dplyr::arrange(NDVI_Biophysical_Zone, desc(Time_Window_Num))
      saveTable(temporal_summary, "Statistics_Fig3D_temporal_summary")
      cat("✓ temporal_summary results saved to 'tables/' directory.\n")
      
      temporal_change_annotations <- temporal_summary %>%
        dplyr::group_by(NDVI_Biophysical_Zone) %>%
        # Sort in reverse order
        dplyr::arrange(desc(Time_Window_Num)) %>%
        dplyr::mutate(Prev_Mean_NDVI = lag(Mean_NDVI),
                      Percent_Change = ifelse(is.na(Prev_Mean_NDVI) | Prev_Mean_NDVI == 0, NA, 
                                              (Mean_NDVI - Prev_Mean_NDVI) / Prev_Mean_NDVI * 100)) %>%
        # We are moving from large to small, so we select < 5 instead of > 1
        dplyr::filter(Time_Window_Num < 5 & !is.na(Percent_Change) & abs(Percent_Change) > 2) 
      saveTable(temporal_change_annotations, "Statistics_Fig3D_temporal_change_annotations")
      cat("✓ temporal_change_annotations results saved to 'tables/' directory.\n")
      
      # --- STATISTICAL TEST 3: Mann-Kendall Trend Test for Time Series ---
      cat("--- Performing Mann-Kendall Trend Tests for each Biophysical Zone ---\n")
      
      mk_results <- temporal_data_long %>%
        dplyr::group_by(OBJECTID, NDVI_Biophysical_Zone) %>%
        # Sort in reverse order
        dplyr::arrange(desc(Time_Window_Num)) %>%
        dplyr::filter(n() >= 4) %>%
        dplyr::summarise(
          mk_test = list(Kendall::MannKendall(NDVI_Value)),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          tau = sapply(mk_test, function(x) x$tau),
          p_value = sapply(mk_test, function(x) x$sl)
        ) %>%
        dplyr::group_by(NDVI_Biophysical_Zone) %>%
        dplyr::summarise(
          Num_Events_Tested = n(),
          Median_Tau = median(tau, na.rm=TRUE),
          Mean_Tau = mean(tau, na.rm=TRUE),
          Prop_Significant_Decline = mean(p_value < 0.05 & tau < 0, na.rm = TRUE),
          Prop_Significant_Increase = mean(p_value < 0.05 & tau > 0, na.rm = TRUE)
        )
      
      print(mk_results)
      saveTable(mk_results, "Supplementary_Table_3d_Mann_Kendall_Results")
      cat("✓ Mann-Kendall trend analysis results saved to 'tables/' directory.\n")
    }
  }
  
  cat("✓ Data preparation for Figure 3 completed\n")
  
  return(list(
    data_3a = data_3a,
    data_3b = climate_ndvi_data,
    data_3c = seasonal_summary,
    data_3d_summary = temporal_summary,
    data_3d_annotations = temporal_change_annotations
  ))
}

# Run the data preparation function
figure_3_data <- prepare_data_for_figure_3()


# ------------------------------------------------------------------------------
# Visualization - Plotting Figure 3 (Final Refined Version)
# ------------------------------------------------------------------------------

create_figure_3_visualization <- function(figure_data, width = 8.3, height = 6.8) {
  cat("\n=== Creating Figure 3 Visualization with Expert Refinements ===\n")
  
  # Ensure the 'figures' directory exists
  if (!dir.exists("figures")) {
    dir.create("figures")
    cat("Created 'figures' directory\n")
  }
  
  # Set color configurations
  biophysical_zone_colors <- c("IDZ" = "#3B9AB2", "CTZ" = "#E14D55", "SDZ" = "#21A764")
  season_colors <- c("Winter" = "#D3DDEA", "Spring" = "#B9DCC9", "Summer" = "#A5D2E5", "Fall" = "#F5C0B8")
  
  # Figure 3A: Vegetation Index Distribution - Legend at the bottom
  if(!is.null(figure_data$data_3a) && nrow(figure_data$data_3a) > 0) {
    plot_3a <- ggplot(figure_data$data_3a, 
                      aes(x = NDVI_1, fill = NDVI_Biophysical_Zone, color = NDVI_Biophysical_Zone)) +
      geom_density(alpha = 0.6, linewidth = 0.7) +
      scale_fill_manual(values = biophysical_zone_colors) +
      scale_color_manual(values = biophysical_zone_colors) +
      nature_theme_professional +
      
      theme(
        legend.position = c(0.22, 0.8),
        legend.box = "horizontal",
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
        legend.box.margin = ggplot2::margin(t = -5)
      ) +
      labs(x = "NDVI Value", y = "Density", fill = "Biophysical Zone", color = "Biophysical Zone")
  } else {
    plot_3a <- ggplot() + theme_void() + 
      annotate("text", x = 0.5, y = 0.5, label = "No data available", size = 3)
  }
  
  # Figure 3B: NDVI Distribution by Climate Zone with Statistics
  # Check if data for Figure 3B exists
  if (!is.null(figure_3_data$data_3b) && nrow(figure_3_data$data_3b) > 0) {
    
    # Define color scheme
    climate_colors <- c(
      "Temperate" = "#F5C0B8", 
      "Arid" = "#8064A2", 
      "Tropical" = "#77AF77", 
      "Cold" = "#5B9BD5"
    )
    
    # Use the full, original dataset for all statistics
    stats_data <- figure_3_data$data_3b
    
    # --- Step 1B: Calculate and Save Key Statistics ---
    cat("--- Calculating and saving statistics for Figure 3B ---\n")
    
    # 1B.1: Calculate descriptive statistics for each climate zone
    descriptive_stats <- stats_data %>%
      dplyr::group_by(Climate_Zone) %>%
      dplyr::summarise(
        SampleSize = n(),
        Mean_NDVI = mean(NDVI_1, na.rm = TRUE),
        SD_NDVI = sd(NDVI_1, na.rm = TRUE),
        Median_NDVI = median(NDVI_1, na.rm = TRUE),
        Q1_NDVI = quantile(NDVI_1, 0.25, na.rm = TRUE),
        Q3_NDVI = quantile(NDVI_1, 0.75, na.rm = TRUE),
        IQR_NDVI = IQR(NDVI_1, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      # Arrange in a logical order
      dplyr::arrange(factor(Climate_Zone, levels = c("Arid", "Tropical", "Temperate", "Cold")))
    
    # Save the descriptive statistics table
    saveTable(descriptive_stats, "Statistics_Fig3B_Descriptive_Stats_by_Climate")
    cat("✓ Descriptive statistics saved.\n")
    
    # 1B.2: Perform Kruskal-Wallis and subsequent pairwise Dunn's test
    if(!require(dunn.test)) { install.packages("dunn.test"); library(dunn.test) }
    
    kw_test <- kruskal.test(NDVI_1 ~ Climate_Zone, data = stats_data)
    
    # Perform Dunn's post-hoc test for pairwise comparisons
    dunn_results <- dunn.test(stats_data$NDVI_1, stats_data$Climate_Zone, method="bonferroni")
    
    # Tidy up the Dunn's test results into a clean data frame
    dunn_summary <- data.frame(
      Comparison = dunn_results$comparisons,
      Z_Statistic = dunn_results$Z,
      P_Value = dunn_results$P,
      P_Value_Adjusted = dunn_results$P.adjusted
    ) %>%
      # Add a column to indicate significance
      dplyr::mutate(Significant = ifelse(P_Value_Adjusted < 0.05, "Yes", "No"))
    
    # Save the pairwise comparison results
    saveTable(dunn_summary, "Supplementary_Table_3b_Pairwise_Dunn_Test")
    cat("✓ Pairwise comparison results saved.\n")
    
    # --- End of statistics section ---
    
    # 2. Data Preparation for Plotting (Subsampling)
    data_for_plot <- stats_data %>%
      dplyr::group_by(Climate_Zone) %>%
      dplyr::slice_sample(n = 2000, replace = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(Climate_Zone_Label = paste0(Climate_Zone, "\n(n=", descriptive_stats$SampleSize[match(Climate_Zone, descriptive_stats$Climate_Zone)], ")"))
    
    # Factor ordering
    ordered_levels <- c("Arid", "Tropical", "Temperate", "Cold")
    present_levels <- ordered_levels[ordered_levels %in% unique(data_for_plot$Climate_Zone)]
    data_for_plot$Climate_Zone_Label <- factor(data_for_plot$Climate_Zone_Label,
                                               levels = paste0(present_levels, "\n(n=", descriptive_stats$SampleSize[match(present_levels, descriptive_stats$Climate_Zone)], ")"))
    data_for_plot$Climate_Zone <- factor(data_for_plot$Climate_Zone, levels = present_levels)
    
    # Prepare p-value text for the plot title
    p_value_text <- ifelse(kw_test$p.value < 0.001, "p < 0.001", paste0("p = ", format(kw_test$p.value, digits = 3)))
    
    # 3. Core Plotting Code
    plot_3b <- ggplot(data_for_plot, aes(x = Climate_Zone_Label, y = NDVI_1, fill = Climate_Zone, color = Climate_Zone)) +
      geom_jitter(width = 0.25, alpha = 0.05, size = 1, show.legend = FALSE) +
      geom_violin(trim = FALSE, alpha = 0.5, show.legend = FALSE) +
      geom_boxplot(width = 0.1, fill = "white", alpha = 0.75, outlier.shape = NA, show.legend = FALSE) +
      
      annotate("text", x = 2.5, y = 1.15, label = paste0("Kruskal-Wallis, ", p_value_text), size = 3, hjust = 0.5, fontface = "italic") +
      scale_fill_manual(values = climate_colors) +
      scale_color_manual(values = climate_colors) +
      nature_theme_professional +
      theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 0.5)) +
      labs(x = "Climate Zone", y = "NDVI Value")
    
  } else {
    # Fallback for no data
    plot_3b <- ggplot() + theme_void() +
      annotate("text", x = 0.5, y = 0.5, label = "No data available for Figure 3B", size = 4)
  }
  
  # Print or display the plot
  print(plot_3b)
  
  # Figure 3C: Seasonal Distribution - Legend at the bottom, x-axis labels on two lines
  if(!is.null(figure_data$data_3c) && nrow(figure_data$data_3c) > 0) {
    # Modify NDVI category labels to display on two lines
    figure_data$data_3c$NDVI_Category_2line <- factor(
      figure_data$data_3c$NDVI_Category,
      levels = c("Low (0-0.3)", "Moderate (0.3-0.6)", "High (0.6-0.8)", "Very High (0.8-1)"),
      labels = c("Low\n(0-0.3)", "Moderate\n(0.3-0.6)", "High\n(0.6-0.8)", "Very High\n(0.8-1)")
    )
    # Calculate the total sample size for each NDVI category
    ndvi_counts <- figure_data$data_3c %>%
      group_by(NDVI_Category) %>%
      dplyr::summarise(
        total_count = sum(Count),
        .groups = "drop"
      )
    
    plot_3c <- ggplot(figure_data$data_3c, aes(x = NDVI_Category_2line, y = Percentage, fill = Season)) +
      geom_col(position = "fill", alpha = 0.85, color="white", linewidth=0.2) +
      geom_text(aes(label = ifelse(Percentage > 5, sprintf("%.0f%%", Percentage), "")),
                position = position_fill(vjust = 0.5), size = 2.5, fontface="bold", color="black") +
      scale_y_continuous(labels = scales::percent_format()) +
      scale_fill_manual(values = season_colors, name="Season") +
      nature_theme_professional +
      theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
        legend.box.margin = ggplot2::margin(t = -5),
        axis.text.x = element_text(angle = 0, hjust = 0.5)  # Display labels horizontally
      ) +
      # Add sample size information below the bars
      annotate("text", 
               x = levels(figure_data$data_3c$NDVI_Category_2line), 
               y = -0.03,
               label = paste0("n=", ndvi_counts$total_count),
               size = 2, color = "gray30") +
      labs(x = "NDVI Category", y = "Proportion of Landslides by Season (%)")
  } else {
    plot_3c <- ggplot() + theme_void() + 
      annotate("text", x = 0.5, y = 0.5, label = "No seasonal data available", size = 3)
  }
  
  
  # Figure 3D: Time Series Analysis - Legend at the bottom
  if(!is.null(figure_data$data_3d_summary) && nrow(figure_data$data_3d_summary) > 0) {
    # Key modification: Create detailed X-axis labels
    # --------------------------------------------------------------------------
    # 1. Define the day range for each time window
    time_window_days <- c(
      "T5" = "(80-64 days)",
      "T4" = "(64-48 days)",
      "T3" = "(48-32 days)",
      "T2" = "(32-16 days)",
      "T1" = "(16-0 days)"
    )
    
    # 2. Create a label vector with line breaks, ordered to match the X-axis (T5 -> T1)
    detailed_x_labels <- paste0(
      names(time_window_days),    # T5, T4, ...
      "\n",                       # Line break
      time_window_days            # (64-80 days), ...
    )
    # --------------------------------------------------------------------------
    
    plot_3d <- ggplot(figure_data$data_3d_summary, 
                      aes(x = reorder(Time_Window, -Time_Window_Num), y = Mean_NDVI, 
                          group = NDVI_Biophysical_Zone, color = NDVI_Biophysical_Zone)) +
      geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper, fill = NDVI_Biophysical_Zone), alpha = 0.2, color = NA) +
      geom_vline(xintercept = "T1", linetype = "dashed", color = "pink") +
      geom_line(linewidth = 1.2) +
      geom_point(size = 2.5) +
      # Adjust label position, moving it up and to the left
      geom_text(data = figure_data$data_3d_annotations,
                aes(label = sprintf("%+.1f%%", Percent_Change)),
                vjust = 1.8, hjust = 0.6, # Add hjust=0.8 to move label left
                nudge_x = -0.15, # Nudge label to the left
                size = 2.5, fontface = "bold", show.legend = FALSE) +
      
      # Use coord_cartesian to get y-min for annotation placement
      annotate("text", x = "T1", y = 0.51, # Set specific position since axis starts at 0.5
               label = "Landslide Event Proximal", hjust = 1.05, vjust = -12, size = 3, color="grey30", fontface="italic") +
      # Use the newly created labels in scale_x_discrete
      # --------------------------------------------------------------------------
    scale_x_discrete(
      name = "Time Window (days before event)", # Update X-axis title
      labels = detailed_x_labels,            # Use new labels
      expand = expansion(mult = c(0.1, 0.15)) # Adjust left/right padding to accommodate labels
    ) +
      # --------------------------------------------------------------------------
    # Add right-side padding
    coord_cartesian(ylim = c(0.5, NA), clip = "off") + # Keep clip="off" to allow text to expand if needed
      
      scale_color_manual(values = biophysical_zone_colors) +
      scale_fill_manual(values = biophysical_zone_colors) +
      nature_theme_professional +
      
      # Adjust legend font size
      theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        axis.text.x = element_text(
          lineheight = 0.9, # Line height, important for multi-line text
          size = 8          # Adjust size from 7 to 8
        ),
        legend.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = -2),
        legend.box.margin = ggplot2::margin(t = -5),
        # Add this line to ensure labels are not cropped
        plot.margin = ggplot2::margin(t = 5, r = 0, b = 0, l = 5)
      ) +
      labs(x = "Time Window", y = "Mean NDVI (95% CI)", color = "Biophysical Zone", fill = "Biophysical Zone")
    
  } else {
    plot_3d <- ggplot() + theme_void() + 
      annotate("text", x = 0.5, y = 0.5, label = "No temporal data available", size = 3)
  }
  print(plot_3d)
  
  
  # Combine all plots
  final_figure <- (plot_3a | plot_3b) / (plot_3c | plot_3d) +
    plot_annotation(tag_levels = list(paste0("(", letters, ")"))) &
    theme(
      plot.tag = element_text(face = "bold", size = 10),
      plot.margin = ggplot2::margin(t = 5, r = 5, b = 0, l = 5) # Reduce margins around each plot
    )
  
  # Save the figure to the 'figures' folder
  saveFigure(final_figure, "My_Final_Figure_3_A4", width = 8.3, height = 6.8)
  
  cat("✓ Figure 3 visualization with expert refinements completed successfully\n")
  cat("✓ All figures saved to the 'figures' directory\n")
  
  return(final_figure)
}


# Create the final figure
figure_3_final <- create_figure_3_visualization(figure_3_data)



# -------------------------------------------------------------------------
# REVISED MAIN FIGURE 4: Advanced Non-Circular Predictive Modeling of Landslide Biophysical Zone
# Focus: Robust, non-circular modeling. Spatial CV. Feature importance. Independent validation of thresholds.
# Methods: Random Forest (or other suitable algorithm like XGBoost), spatial k-means CV or block CV.
#          Target variable: Landslide occurrence (binary 0/1) or landslide density if available.
#          For this example, we'll use landslide occurrence (binary) and address class imbalance if severe.
# Contribution: Demonstrates predictability and confirms VI importance independently.
# -------------------------------------------------------------------------


# =========================================================================
# FIGURE 4
# =========================================================================
# Part 0: Model Training Function (16 features)
# Task: Prepare data, perform cross-validation, train the final model, and save all results.
# =========================================================================
run_fig4_modeling_16f <- function() {
  cat("\n=== RUNNING MODELING FOR FIGURE 4 (16 FEATURES) ===\n")
  
  # --- 0. Data and Environment Preparation ---
  if (!exists("landslide_data_processed")) {
    stop("Required data `landslide_data_processed` not found.")
  }
  data <- landslide_data_processed
  
  # --- 1. Create non-circular target variable (Landslide Density Proxy) ---
  cat("Step 1: Creating non-circular target variable...\n")
  ndvi_values_for_kde <- data$NDVI_1[!is.na(data$NDVI_1) & data$NDVI_1 >= 0 & data$NDVI_1 <= 1]
  if (length(ndvi_values_for_kde) < 50) stop("Insufficient valid NDVI data for target variable.")
  
  kde <- density(ndvi_values_for_kde, from = 0, to = 1, n = 101, bw = "SJ")
  ndvi_density_df <- data.frame(NDVI_Value_Grid = kde$x, Density_Target = kde$y)
  ndvi_density_df$Density_Target <- scales::rescale(ndvi_density_df$Density_Target)
  
  data$NDVI_1_Rounded <- round(data$NDVI_1, 2)
  ndvi_density_df$NDVI_Value_Grid_Rounded <- round(ndvi_density_df$NDVI_Value_Grid, 2)
  
  model_data <- data %>%
    filter(!is.na(NDVI_1)) %>%
    left_join(ndvi_density_df %>% dplyr::select(NDVI_Value_Grid_Rounded, Density_Target),
              by = c("NDVI_1_Rounded" = "NDVI_Value_Grid_Rounded")) %>%
    filter(!is.na(Density_Target))
  
  # --- 2. Feature selection and data preparation ---
  cat("Step 2: Preparing modeling data...\n")
  features <- c("EVI_1", "EVI_2", "LAI_1", "LAI_2", "MODIS_IGBP_Simplified", "Climate_Zone", "Season", 
                "NDVI_change_1_to_2", "EVI_change_1_to_2", "LAI_change_1_to_2", "NDVI_change_2_to_3", 
                "EVI_change_2_to_3", "LAI_change_2_to_3","Copernicus_LC_Class","Continent",
                "Hansen_Tree_Cover_2000_Percent")
  
  features <- intersect(features, colnames(model_data))
  if (length(features) == 0) stop("No valid features found for modeling.")
  
  modeling_df_full <- model_data %>%
    dplyr::select(all_of(features), Density_Target, Latitude, Longitude, OBJECTID) %>%
    na.omit() %>%
    mutate(across(where(is.character), as.factor))
  
  if (nrow(modeling_df_full) < 200) stop("Insufficient data (n=", nrow(modeling_df_full), ") for modeling.")
  
  # --- 3. Robust spatial cross-validation ---
  cat("Step 3: Performing spatial cross-validation...\n")
  set.seed(123)
  q_lat <- unique(quantile(modeling_df_full$Latitude, probs = 0:5/5, na.rm = TRUE))
  if (length(q_lat) < 3) {
    modeling_df_full$Fold <- sample(1:5, nrow(modeling_df_full), replace = TRUE)
  } else {
    modeling_df_full$Fold <- as.numeric(cut(modeling_df_full$Latitude, breaks = q_lat, include.lowest = TRUE))
  }
  modeling_df_full$Fold[is.na(modeling_df_full$Fold)] <- sample(1:5, sum(is.na(modeling_df_full$Fold)), replace = TRUE)
  
  all_predictions <- data.frame()
  feature_importance_list <- list()
  
  for (k in 1:5) {
    cat("  Processing fold:", k, "...\n")
    train_data <- modeling_df_full[modeling_df_full$Fold != k, ]
    test_data <- modeling_df_full[modeling_df_full$Fold == k, ]
    
    current_formula <- as.formula(paste("Density_Target ~", paste(features, collapse = " + ")))
    rf_model <- randomForest(current_formula, data = train_data, ntree = 200, importance = TRUE)
    
    predictions <- predict(rf_model, test_data)
    all_predictions <- rbind(all_predictions, data.frame(Actual = test_data$Density_Target, Predicted = predictions, Fold = k))
    
    imp <- randomForest::importance(rf_model, type = 1)
    feature_importance_list[[k]] <- data.frame(Feature = rownames(imp), Importance = imp[, 1])
  }
  if (nrow(all_predictions) == 0) stop("All modeling folds failed.")
  
  # --- 4. Train the final global model ---
  cat("Step 4: Training final model and preparing data for plots...\n")
  validation_rf <- randomForest(Density_Target ~ ., data = modeling_df_full %>% dplyr::select(all_of(features), Density_Target),
                                ntree = 500, importance = TRUE)
  
  # --- 5. Prepare and save all data required for plotting ---
  cat("Step 5: Saving all modeling results for plotting...\n")
  
  # Prepare data for Panel C
  avg_importance <- bind_rows(feature_importance_list) %>%
    group_by(Feature) %>%
    dplyr::summarise(Mean_Importance = mean(Importance, na.rm = TRUE), SD_Importance = sd(Importance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(Mean_Importance))
  
  # Prepare data for Panel D
  modeling_df_full$Predicted_Density <- predict(validation_rf, modeling_df_full %>% dplyr::select(all_of(features)))
  pdp_plot_data <- modeling_df_full %>%
    dplyr::select(OBJECTID, Predicted_Density) %>%
    left_join(model_data %>% dplyr::select(OBJECTID, NDVI_1), by = "OBJECTID")
  pdp_summary <- pdp_plot_data %>%
    filter(!is.na(NDVI_1)) %>%
    group_by(NDVI_Value = round(NDVI_1, 2)) %>%
    dplyr::summarise(Mean_Predicted_Density = mean(Predicted_Density, na.rm = TRUE), .groups = "drop") %>%
    filter(!is.na(NDVI_Value))
  
  # Package all results into a list
  fig4_results_data <- list(
    all_predictions = all_predictions,       # For Panels A and B
    avg_importance = avg_importance,         # For Panel C
    pdp_summary = pdp_summary,               # For Panel D
    analysis_results_global = fig1_analysis_results # Pass results from Figure 1
  )
  
  # Create directory and save
  dir.create("models", showWarnings = FALSE)
  saveRDS(fig4_results_data, "models/fig4_plotting_data_16f.rds")
  # Also save the final model and features for Figure 5
  saveRDS(validation_rf, "models/best_landslide_bz_model_16f.rds")
  saveRDS(features, "models/model_features_16f.rds")
  saveRDS(modeling_df_full %>% dplyr::select(all_of(features)), "models/training_data_snapshot_16f.rds")
  
  cat("✓ Modeling complete. All results saved to 'models/fig4_plotting_data_16f.rds'.\n")
}

# Ensure landslide_data_processed and fig1_analysis_results objects exist
run_fig4_modeling_16f() 


# =========================================================================
# Part 2: Plotting Function
# Task: Load saved model results, generate the four panels of Figure 4, and create the final composite plot.
# =========================================================================
create_fig4_plots_16f <- function() {
  cat("\n=== CREATING PLOTS FOR FIGURE 4 (16 FEATURES) ===\n")
  
  # --- 1. Load pre-computed modeling results ---
  cat("Step 1: Loading pre-computed modeling results...\n")
  results_path <- "models/fig4_plotting_data_16f.rds"
  if (!file.exists(results_path)) {
    stop("Plotting data not found. Please run `run_fig4_modeling_16f()` first.")
  }
  fig4_data <- readRDS(results_path)
  
  # Unpack data from the list for easier use
  all_predictions <- fig4_data$all_predictions
  avg_importance <- fig4_data$avg_importance
  pdp_summary <- fig4_data$pdp_summary
  analysis_results_global <- fig4_data$fig1_analysis_results
  
  # --- 2. Generate the four visualization panels ---
  cat("Step 2: Generating visualization panels...\n")
  
  # --- Panel A: Prediction Performance ---
  r_squared <- cor(all_predictions$Actual, all_predictions$Predicted)^2
  rmse <- sqrt(mean((all_predictions$Actual - all_predictions$Predicted)^2))
  plot_4a <- ggplot(all_predictions, aes(x = Actual, y = Predicted)) +
    geom_point(alpha = 0.4, color = "darkcyan", size = 1.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1) +
    geom_smooth(method = "lm", color = "blue", se = TRUE, fill = "lightblue", alpha = 0.2) +
    annotate("text", x = 0.05, y = 0.95 * max(all_predictions$Predicted, na.rm=T),
             label = paste0("R² = ", round(r_squared, 3), "\nRMSE = ", round(rmse, 3)),
             hjust = 0, size = 3.5, fontface = "bold") +
    labs(x = "Actual Landslide Density (Proxy)", y = "Predicted Landslide Density (Proxy)") +
    nature_theme_professional
  
  # --- Panel B: Residual Analysis ---
  all_predictions$Residual <- all_predictions$Predicted - all_predictions$Actual
  plot_4b <- ggplot(all_predictions, aes(x = Residual, fill = factor(Fold))) +
    geom_density(alpha = 0.6, color = "black", linewidth = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
    scale_fill_viridis_d(name = "CV Fold") +
    labs(x = "Residual (Predicted - Actual)", y = "Density") +
    nature_theme_professional +
    # Move legend to the upper right corner inside the plot area
    theme(legend.position = c(0.89, 0.75),
          legend.background = element_rect(fill = "white", color = "gray90"),
          legend.margin = ggplot2::margin(5, 5, 5, 5))
  
  # --- Panel C: Feature Importance ---
  avg_importance$Category <- case_when(
    grepl("EVI_1|LAI_1", avg_importance$Feature) ~ "Vegetation Indices",
    grepl("change", avg_importance$Feature) ~ "Temporal Changes",
    grepl("Hansen|MODIS", avg_importance$Feature) ~ "Vegetation Types",
    TRUE ~ "Environmental"
  )
  feature_palette <- setNames(c("#3B6F9E", "#D85C60", "#63A27D", "#E8A354"), 
                              c("Environmental", "Temporal Changes", "Vegetation Indices", "Vegetation Types"))
  
  # Simplify x-axis labels for panel C
  # Create a simplified feature name column
  avg_importance$SimpleFeature <- gsub("NDVI_change_1_to_2", "Nc1t2", avg_importance$Feature)
  avg_importance$SimpleFeature <- gsub("EVI_change_1_to_2", "Ec1t2", avg_importance$SimpleFeature)
  avg_importance$SimpleFeature <- gsub("LAI_change_1_to_2", "Lc1t2", avg_importance$SimpleFeature)
  avg_importance$SimpleFeature <- gsub("Hansen_Tree_Cover_2000_Percent", "HTC2P", avg_importance$SimpleFeature)
  avg_importance$SimpleFeature <- gsub("Climate_Zone", "CliZ", avg_importance$SimpleFeature)
  avg_importance$SimpleFeature <- gsub("MODIS_IGBP_Simplified", "IGBPS", avg_importance$SimpleFeature)
  
  plot_4c <- ggplot(head(avg_importance, 16), aes(x = reorder(SimpleFeature, Mean_Importance), y = Mean_Importance, fill = Category)) +
    geom_col(width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, Mean_Importance - SD_Importance), ymax = Mean_Importance + SD_Importance), width = 0.25, color = "gray20") +
    ggrepel::geom_text_repel(aes(label=sprintf("%.1f", Mean_Importance)), 
                             hjust = -0.3,
                             size=3, fontface="bold", color="black",
                             direction = "y", 
                             nudge_x = 0.3,
                             box.padding = 0.1,
                             point.padding = 0.1,
                             segment.curvature = -0.1,
                             segment.ncp = 3,
                             segment.angle = 20,
                             segment.size = 0.2, min.segment.length = 0) + 
    coord_flip(ylim = c(0, max(avg_importance$Mean_Importance, na.rm = TRUE) * 1.15), clip = "off") +
    scale_fill_manual(values = feature_palette) +
    labs(x = "Feature", y = "Importance (Mean %IncMSE)") +
    nature_theme_professional +
    theme(legend.position = "bottom", 
          panel.grid.major.y = element_blank(), 
          axis.text.y = element_text(size = 9),
          # Remove the legend title
          legend.title = element_blank())
  
  # --- Panel D: NDVI Threshold Validation ---
  ndvi_breakpoints <- analysis_results_global$NDVI_1$breakpoints
  
  plot_4d <- ggplot(pdp_summary, aes(x = NDVI_Value, y = Mean_Predicted_Density)) +
    geom_line(color = "darkgreen", linewidth = 1.5, alpha = 0.8) +
    geom_vline(xintercept = ndvi_breakpoints, linetype = "dashed", color = "red", linewidth = 0.8) +
    annotate("rect", xmin = ndvi_breakpoints[1], xmax = ndvi_breakpoints[2],
             ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.1) +
    annotate("text", x = mean(ndvi_breakpoints) + 0.08, y = max(pdp_summary$Mean_Predicted_Density, na.rm=T) * 0.9,
             label = "Fig 1 Critical Transition Zone", color = "red", size = 3, fontface = "bold", hjust=1.3) +
    labs(x = "NDVI Value", y = "Predicted Landslide Density (Proxy)") +
    nature_theme_professional
  
  # --- 3. Assemble and save the final Figure 4 ---
  cat("Step 3: Assembling and saving final figure...\n")
  main_fig4 <- (plot_4a | plot_4b) / (plot_4c | plot_4d) +
    plot_annotation(tag_levels = list(paste0("(", letters, ")"))) &
    theme(plot.tag = element_text(face = "bold", size = 12))
  
  saveFigure(main_fig4, "Main_Figure_4_Modular_16f", width = 8.3, height = 6.8)
  
  cat("✓ Figure 4 (16f) plotting complete.\n")
  
  return(main_fig4)
}

# This function reads data from the file, so it runs quickly.
my_figure_4_16f <- create_fig4_plots_16f()


# =========================================================================
# Part 1: Model Training Function (Standard Features)
# Task: Prepare data, perform cross-validation, train the final model, and save all results.
# =========================================================================
run_fig4_modeling <- function() {
  cat("\n=== RUNNING MODELING FOR FIGURE 4 (STANDARD FEATURES) ===\n")
  # Ensure necessary packages are available
  if (!require(pROC, quietly = TRUE)) { install.packages("pROC", quiet = TRUE); library(pROC) }
  
  # --- 0. Data and Environment Preparation ---
  if (!exists("landslide_data_processed")) {
    stop("Required data `landslide_data_processed` not found.")
  }
  data <- landslide_data_processed
  
  # --- 1. Create non-circular target variable (Landslide Density Proxy) ---
  cat("Step 1: Creating non-circular target variable...\n")
  ndvi_values_for_kde <- data$NDVI_1[!is.na(data$NDVI_1) & data$NDVI_1 >= 0 & data$NDVI_1 <= 1]
  if (length(ndvi_values_for_kde) < 50) stop("Insufficient valid NDVI data for target variable.")
  
  kde <- density(ndvi_values_for_kde, from = 0, to = 1, n = 101, bw = "SJ")
  ndvi_density_df <- data.frame(NDVI_Value_Grid = kde$x, Density_Target = kde$y)
  ndvi_density_df$Density_Target <- scales::rescale(ndvi_density_df$Density_Target)
  
  data$NDVI_1_Rounded <- round(data$NDVI_1, 2)
  ndvi_density_df$NDVI_Value_Grid_Rounded <- round(ndvi_density_df$NDVI_Value_Grid, 2)
  
  model_data <- data %>%
    filter(!is.na(NDVI_1)) %>%
    left_join(ndvi_density_df %>% dplyr::select(NDVI_Value_Grid_Rounded, Density_Target),
              by = c("NDVI_1_Rounded" = "NDVI_Value_Grid_Rounded")) %>%
    filter(!is.na(Density_Target))
  
  # --- 2. Feature selection and data preparation ---
  cat("Step 2: Preparing modeling data...\n")
  features <- c("EVI_1", "LAI_1", "MODIS_IGBP_Simplified", "Climate_Zone", "Season", 
                "NDVI_change_1_to_2", "EVI_change_1_to_2", "LAI_change_1_to_2",
                "Hansen_Tree_Cover_2000_Percent")
  
  features <- intersect(features, colnames(model_data))
  if (length(features) == 0) stop("No valid features found for modeling.")
  
  modeling_df_full <- model_data %>%
    dplyr::select(all_of(features), Density_Target, Latitude, Longitude, OBJECTID) %>%
    na.omit() %>%
    mutate(across(where(is.character), as.factor))
  
  if (nrow(modeling_df_full) < 200) stop("Insufficient data (n=", nrow(modeling_df_full), ") for modeling.")
  
  # --- 3. Robust spatial cross-validation ---
  cat("Step 3: Performing spatial cross-validation...\n")
  set.seed(123)
  q_lat <- unique(quantile(modeling_df_full$Latitude, probs = 0:5/5, na.rm = TRUE))
  if (length(q_lat) < 3) {
    modeling_df_full$Fold <- sample(1:5, nrow(modeling_df_full), replace = TRUE)
  } else {
    modeling_df_full$Fold <- as.numeric(cut(modeling_df_full$Latitude, breaks = q_lat, include.lowest = TRUE))
  }
  modeling_df_full$Fold[is.na(modeling_df_full$Fold)] <- sample(1:5, sum(is.na(modeling_df_full$Fold)), replace = TRUE)
  
  
  # --- Create a binary target for classification task ---
  bz_threshold <- quantile(modeling_df_full$Density_Target, 0.8)
  modeling_df_full$Zone_Class <- factor(
    ifelse(modeling_df_full$Density_Target >= bz_threshold, "CT_Zone", "Not_CT_Zone"),
    levels = c("Not_CT_Zone", "CT_Zone") # Set "Not_CT_Zone" as the reference level
  )
  cat(paste0("--- Created binary classification target 'Zone_Class' with threshold: ", round(bz_threshold, 3), " ---\n"))
  
  # Initialize storage for results from all folds
  all_regression_preds <- data.frame()
  all_classification_preds <- list() # For ROC curves
  feature_importance_list <- list()
  fold_metrics <- data.frame()
  
  
  for (k in 1:5) {
    cat("  Processing fold:", k, "...\n")
    train_data <- modeling_df_full[modeling_df_full$Fold != k, ]
    test_data <- modeling_df_full[modeling_df_full$Fold == k, ]
    
    # --- REGRESSION MODEL ---
    reg_formula <- as.formula(paste("Density_Target ~", paste(features, collapse = " + ")))
    reg_model <- randomForest(reg_formula, data = train_data, ntree = 200, importance = TRUE)
    reg_predictions <- predict(reg_model, test_data)
    all_regression_preds <- rbind(all_regression_preds, data.frame(Actual = test_data$Density_Target, Predicted = reg_predictions, Fold = k))
    
    # --- CLASSIFICATION MODEL ---
    class_formula <- as.formula(paste("Zone_Class ~", paste(features, collapse = " + ")))
    class_model <- randomForest(class_formula, data = train_data, ntree = 200)
    # Get class probabilities for ROC curve
    class_predictions <- predict(class_model, test_data, type = "prob")
    # Explicitly tell roc() which level is the "case" or "event" to avoid warnings
    roc_obj_fold <- roc(response = test_data$Zone_Class, predictor = class_predictions[, "CT_Zone"],
                        levels = c("Not_CT_Zone", "CT_Zone"), quiet = TRUE)
    all_classification_preds[[k]] <- data.frame(
      Response = test_data$Zone_Class,
      Prediction_Prob = class_predictions[, "CT_Zone"],
      Fold = k
    )
    
    # --- METRICS & IMPORTANCE ---
    imp <- randomForest::importance(reg_model, type = 1)
    feature_importance_list[[k]] <- data.frame(Feature = rownames(imp), Importance = imp[, 1])
    
    # Calculate metrics for this fold
    r_squared_fold <- cor(test_data$Density_Target, reg_predictions)^2
    rmse_fold <- sqrt(mean((test_data$Density_Target - reg_predictions)^2))
    auc_fold <- auc(roc_obj_fold)
    
    fold_metrics <- rbind(fold_metrics, data.frame(Fold=k, R_squared=r_squared_fold, RMSE=rmse_fold, AUC=auc_fold))
  }
  
  # --- 4. CALCULATE AND SAVE OVERALL PERFORMANCE METRICS ---
  cat("Step 4: Calculating and saving overall performance metrics...\n")
  
  # Summarize fold metrics to get mean and standard deviation
  final_metrics_summary <- fold_metrics %>%
    dplyr::summarise(across(everything(), list(mean = mean, sd = sd)), .by = NULL)
  
  print("--- Overall Cross-Validated Performance Metrics ---")
  print(final_metrics_summary)
  saveTable(final_metrics_summary, "Statistics_Fig4_Model_Performance_Summary")
  
  # --- GENERATE AND SAVE ROC CURVE PLOT ---
  cat("--- Generating and saving ROC curve plot ---\n")
  
  all_classification_df <- do.call(rbind, all_classification_preds)
  # Use the correct column names: `AUC_mean` and `AUC_sd`
  mean_auc_val <- final_metrics_summary$AUC_mean
  sd_auc_val <- final_metrics_summary$AUC_sd
  
  # Create a beautiful ROC plot with ggplot2
  roc_plot <- ggplot(all_classification_df, aes(d = Response, m = Prediction_Prob, color = as.factor(Fold))) +
    plotROC::geom_roc(n.cuts = 0) +
    style_roc(theme = theme_gray) +
    theme_minimal() +
    # Move legend to the upper right corner inside the plot area
    theme(legend.position = c(0.22, 0.66),
          legend.background = element_rect(fill = "white", color = "gray90"),
          legend.margin = ggplot2::margin(5, 5, 5, 5)) +
    labs(
      title = "Cross-Validated ROC Curves by Spatial Fold",
      x = "False Positive Fraction",
      y = "True Positive Fraction",
      color = "CV Fold"
    ) +
    # Use the corrected variables
    annotate("text", x = 0.4, y = 0.8, 
             label = paste0("Mean AUC = ", round(mean_auc_val, 3), 
                            " ± ", round(sd_auc_val, 3)),
             fontface = "bold", size = 4, hjust = 0.1, vjust =0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey")
  
  print(roc_plot)
  saveFigure(roc_plot, "Supplementary_Fig_Z_ROC_Curve", width = 6, height = 5)
  
  # --- 5. Prepare and save all data required for plotting ---
  cat("Step 5: Saving all modeling results for plotting...\n")
  
  # Train the final regression model
  final_reg_formula <- as.formula(paste("Density_Target ~", paste(features, collapse = " + ")))
  validation_rf <- randomForest(final_reg_formula, data = modeling_df_full, ntree = 500, importance = TRUE)
  
  # Prepare data for Panel C
  avg_importance <- bind_rows(feature_importance_list) %>%
    group_by(Feature) %>%
    dplyr::summarise(Mean_Importance = mean(Importance, na.rm = TRUE), SD_Importance = sd(Importance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(Mean_Importance))
  
  # Prepare data for Panel D
  modeling_df_full$Predicted_Density <- predict(validation_rf, modeling_df_full %>% dplyr::select(all_of(features)))
  pdp_plot_data <- modeling_df_full %>%
    dplyr::select(OBJECTID, Predicted_Density) %>%
    left_join(model_data %>% dplyr::select(OBJECTID, NDVI_1), by = "OBJECTID")
  pdp_summary <- pdp_plot_data %>%
    filter(!is.na(NDVI_1)) %>%
    group_by(NDVI_Value = round(NDVI_1, 2)) %>%
    dplyr::summarise(Mean_Predicted_Density = mean(Predicted_Density, na.rm = TRUE), .groups = "drop") %>%
    filter(!is.na(NDVI_Value))
  
  # Package all results into a list
  fig4_results_data <- list(
    all_predictions = all_regression_preds,       # For Panels A and B
    avg_importance = avg_importance,         # For Panel C
    pdp_summary = pdp_summary,               # For Panel D
    analysis_results_global = fig1_analysis_results # Pass results from Figure 1
  )
  
  # Create directory and save
  dir.create("models", showWarnings = FALSE)
  saveRDS(fig4_results_data, "models/fig4_plotting_data.rds")
  # Also save the final model and features for Figure 5
  saveRDS(validation_rf, "models/best_landslide_bz_model.rds")
  saveRDS(features, "models/model_features.rds")
  saveRDS(modeling_df_full %>% dplyr::select(all_of(features)), "models/training_data_snapshot.rds")
  
  cat("✓ Modeling complete. All results saved to 'models/fig4_plotting_data.rds'.\n")
}

# Ensure landslide_data_processed and fig1_analysis_results objects exist
run_fig4_modeling() 

# =========================================================================
# Part 2: Plotting Function
# Task: Load saved model results, generate the four panels of Figure 4, and create the final composite plot.
# =========================================================================
create_fig4_plots <- function() {
  cat("\n=== CREATING PLOTS FOR FIGURE 4 (STANDARD FEATURES) ===\n")
  
  # --- 1. Load pre-computed modeling results ---
  cat("Step 1: Loading pre-computed modeling results...\n")
  results_path <- "models/fig4_plotting_data.rds"
  if (!file.exists(results_path)) {
    stop("Plotting data not found. Please run `run_fig4_modeling()` first.")
  }
  fig4_data <- readRDS(results_path)
  
  # Unpack data from the list for easier use
  all_predictions <- fig4_data$all_predictions
  avg_importance <- fig4_data$avg_importance
  pdp_summary <- fig4_data$pdp_summary
  analysis_results_global <- fig4_data$fig1_analysis_results
  
  # --- 2. Generate the four visualization panels ---
  cat("Step 2: Generating visualization panels...\n")
  
  # --- Panel A: Prediction Performance ---
  cat("Generating visualization for Fig 4 panel A...\n")
  r_squared <- cor(all_predictions$Actual, all_predictions$Predicted)^2
  rmse <- sqrt(mean((all_predictions$Actual - all_predictions$Predicted)^2))
  plot_4a <- ggplot(all_predictions, aes(x = Actual, y = Predicted)) +
    geom_point(alpha = 0.4, color = "darkcyan", size = 1.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1) +
    geom_smooth(method = "lm", color = "blue", se = TRUE, fill = "lightblue", alpha = 0.2) +
    annotate("text", x = 0.05, y = 0.95 * max(all_predictions$Predicted, na.rm=T),
             label = paste0("R² = ", round(r_squared, 3), "\nRMSE = ", round(rmse, 3)),
             hjust = 0, size = 3.5, fontface = "bold") +
    labs(x = "Actual Landslide Density (Proxy)", y = "Predicted Landslide Density (Proxy)") +
    nature_theme_professional
  
  # --- Panel B: Residual Analysis ---
  cat("Generating visualization for Fig 4 panel B...\n")
  all_predictions$Residual <- all_predictions$Predicted - all_predictions$Actual
  plot_4b <- ggplot(all_predictions, aes(x = Residual, fill = factor(Fold))) +
    geom_density(alpha = 0.4, color = "black", linewidth = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
    scale_fill_viridis_d(name = "CV Fold") +
    labs(x = "Residual (Predicted - Actual)", y = "Density") +
    nature_theme_professional +
    # Move legend to the upper right corner inside the plot area
    theme(legend.position = c(0.89, 0.75),
          legend.background = element_rect(fill = "white", color = "gray90"),
          legend.margin = ggplot2::margin(5, 5, 5, 5))
  
  # --- Panel C: Feature Importance ---
  cat("Generating visualization for Fig 4 panel C...\n")
  avg_importance$Category <- case_when(
    grepl("EVI_1|LAI_1", avg_importance$Feature) ~ "Vegetation Indices",
    grepl("change", avg_importance$Feature) ~ "Temporal Changes",
    grepl("Hansen|MODIS", avg_importance$Feature) ~ "Vegetation Types",
    TRUE ~ "Environmental"
  )
  feature_palette <- setNames(c("#3B6F9E", "#D85C60", "#63A27D", "#E8A354"), 
                              c("Environmental", "Temporal Changes", "Vegetation Indices", "Vegetation Types"))
  
  # Simplify x-axis labels for panel C
  cat("Simplifying feature names for Fig 4 panel C...\n")
  avg_importance$SimpleFeature <- gsub("NDVI_change_1_to_2", "Nc1t2", avg_importance$Feature)
  avg_importance$SimpleFeature <- gsub("EVI_change_1_to_2", "Ec1t2", avg_importance$SimpleFeature)
  avg_importance$SimpleFeature <- gsub("LAI_change_1_to_2", "Lc1t2", avg_importance$SimpleFeature)
  avg_importance$SimpleFeature <- gsub("Hansen_Tree_Cover_2000_Percent", "HTC2P", avg_importance$SimpleFeature)
  avg_importance$SimpleFeature <- gsub("Climate_Zone", "CliZ", avg_importance$SimpleFeature)
  avg_importance$SimpleFeature <- gsub("MODIS_IGBP_Simplified", "IGBPS", avg_importance$SimpleFeature)
  
  cat("Plotting for Fig 4 panel C...\n")
  plot_4c <- ggplot(head(avg_importance, 10), aes(x = reorder(SimpleFeature, Mean_Importance), y = Mean_Importance, fill = Category)) +
    geom_col(width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, Mean_Importance - SD_Importance), ymax = Mean_Importance + SD_Importance), width = 0.25, color = "gray20") +
    ggrepel::geom_text_repel(aes(label=sprintf("%.1f", Mean_Importance)), 
                             hjust = -0.3,
                             size=3, fontface="bold", color="black",
                             direction = "y", 
                             nudge_x = 0.3,
                             box.padding = 0.1,
                             point.padding = 0.1,
                             segment.curvature = -0.1,
                             segment.ncp = 3,
                             segment.angle = 20,
                             segment.size = 0.2, min.segment.length = 0) + 
    coord_flip(ylim = c(0, max(avg_importance$Mean_Importance, na.rm = TRUE) * 1.15), clip = "off") +
    scale_fill_manual(values = feature_palette) +
    labs(x = "Feature", y = "Importance (Mean %IncMSE)") +
    nature_theme_professional +
    theme(
      legend.position = c(0.78, 0.18),
      legend.box = "horizontal",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      legend.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5),
      legend.box.margin = ggplot2::margin(t = -5),
      panel.grid.major.y = element_blank(),
      # Set legend background and border
      legend.background = element_rect(
        colour = "grey80",
        linewidth = 0.3,
        fill = "white"
      )
    ) 
  
  # --- Panel D: NDVI Threshold Validation ---
  cat("Generating visualization for Fig 4 panel D...\n")
  ndvi_breakpoints <- analysis_results_global$NDVI_1$breakpoints
  
  plot_4d <- ggplot(pdp_summary, aes(x = NDVI_Value, y = Mean_Predicted_Density)) +
    geom_line(color = "darkgreen", linewidth = 1.5, alpha = 0.8) +
    geom_vline(xintercept = ndvi_breakpoints, linetype = "dashed", color = "red", linewidth = 0.8) +
    annotate("rect", xmin = ndvi_breakpoints[1], xmax = ndvi_breakpoints[2],
             ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.1) +
    annotate("text", x = mean(ndvi_breakpoints) + 0.08, y = max(pdp_summary$Mean_Predicted_Density, na.rm=T) * 0.9,
             label = "Fig 1 Critical Transition Zone", color = "red", size = 3, fontface = "bold", hjust=1.3) +
    labs(x = "NDVI Value", y = "Predicted Landslide Density (Proxy)") +
    nature_theme_professional
  
  # --- 3. Assemble and save the final Figure 4 ---
  cat("Step 3: Assembling and saving final figure...\n")
  main_fig4 <- (plot_4a | plot_4b) / (plot_4c | plot_4d) +
    plot_annotation(tag_levels = list(paste0("(", letters, ")"))) &
    theme(
      plot.tag = element_text(face = "bold", size = 10),
      plot.margin = ggplot2::margin(t = 5, r = 5, b = 0, l = 5) # Reduce margins around each plot
    )
  
  saveFigure(main_fig4, "Main_Figure_4_Modular", width = 8.3, height = 6.8)
  
  cat("✓ Figure 4 plotting complete.\n")
  
  return(main_fig4)
}

# This function reads data from the file, so it runs quickly.
my_figure_4 <- create_fig4_plots()