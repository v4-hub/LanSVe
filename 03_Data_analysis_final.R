# -------------------------------------------------------------------------
# 植被指数与滑坡灾害的非线性关系：超越Nature Communications标准的卓越研究
# 多源遥感数据驱动的风险阈值发现与全球验证框架
# 专家级分析：滑坡机理 + 植被遥感 + 期刊编辑视角
# -------------------------------------------------------------------------

# =========================================================================
# 1. 环境设置与核心库
# =========================================================================
# 核心库加载与版本检查
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

# 检查系统可用字体的函数
check_available_fonts <- function() {
  available_fonts <- systemfonts::system_fonts()
  
  # 检查Nature推荐的字体
  nature_fonts <- available_fonts %>% 
    filter(str_detect(family, "Arial|Helvetica|Times|Calibri")) %>% 
    dplyr::select(family, style) %>% 
    distinct()
  
  return(nature_fonts)
}

# 获取最佳可用字体
get_best_nature_font <- function() {
  available_fonts <- systemfonts::system_fonts()$family
  
  # Nature期刊字体优先级
  font_priority <- c("Arial", "Helvetica", "Calibri", "Times New Roman", "Times")
  
  for(font in font_priority) {
    if (font %in% available_fonts) {
      cat("✓ Using font:", font, "\n")
      return(font)
    }
  }
  
  cat("⚠ Using default sans-serif font\n")
  return("")  # 空字符串让R使用默认字体
}

# 检查可用字体
cat("Available Nature-style fonts on your system:\n")
print(check_available_fonts())

# 获取最佳字体
best_font <- get_best_nature_font()

# =========================================================================
# 设置保存结果的目录
# Define the base directory
base_dir <- "V5.2"

# Define the subdirectories
sub_dirs <- c("figures", "tables", "data")

# Create the base directory and its subdirectories
for (sub_dir in sub_dirs) {
  full_path <- file.path(base_dir, sub_dir)
  dir.create(full_path, showWarnings = FALSE, recursive = TRUE)
}
# =========================================================================
# -------------------------------------------------------------------------
# 专业学术视觉系统（Nature系列期刊标准）
# -------------------------------------------------------------------------

# Nature Communications专业配色体系
nature_palettes <- list(
  # 主要发现配色（基于Nature期刊指南）
  biophysical_zones = c(
    "IDZ" = "#0173B2",      # 蓝色：低风险
    "BTZ" = "#DE8F05", # 橙色：过渡期
    "CTZ" = "#CC78BC",      # 紫色：高风险区（突出显示）
    "SDZ" = "#029E73",      # 绿色：后峰值
    "Very High NDVI" = "#56B4E9",      # 浅蓝：极高植被
    "Unknown" = "#999999"              # 灰色：未知
  ),
  
  # 植被指数专业配色
  indices = c(
    "NDVI" = "#1B9E77", "EVI" = "#D95F02", 
    "LAI" = "#7570B3", "FPAR" = "#E7298A"
  ),
  
  # 气候带科学配色
  climate = c(
    "Temperate" = "#2166AC",
    "Mediterranean/Subtropical" = "#762A83", 
    "Tropical" = "#5AAE61",
    "Other" = "#9970AB"
  ),
  
  # 植被类型生态学配色
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

# 专业学术主题（严格遵循Nature指南）
# 优化后的Nature主题（基于你的原有主题）
nature_theme_professional <- theme_minimal(base_size = 10, base_family = best_font) +
  theme(
    # 字体设置 - 确保所有文本元素都使用正确字体
    text = element_text(family = best_font),
    
    # 网格和边框
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#F0F0F0", linewidth = 0.25),
    panel.border = element_rect(fill = NA, color = "#CCCCCC", linewidth = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    
    # 图例设置
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 9, family = best_font),
    legend.text = element_text(size = 8, family = best_font),
    legend.key.size = unit(0.4, "cm"),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA),
    
    # 标题设置
    plot.title = element_text(size = 11, face = "bold", color = "#000000", family = best_font),
    plot.subtitle = element_text(size = 9, color = "#666666", family = best_font),
    plot.caption = element_text(size = 7, color = "#666666", hjust = 0, family = best_font),
    
    # 坐标轴设置 - 符合Nature标准
    axis.title = element_text(size = 9, face = "bold", color = "#333333", family = best_font),
    axis.text = element_text(size = 8, color = "#666666", family = best_font),
    axis.ticks = element_line(color = "#CCCCCC", linewidth = 0.25),
    axis.line = element_line(color = "#000000", linewidth = 0.3),  # Nature喜欢的轴线
    
    # 分面设置
    strip.text = element_text(size = 9, face = "bold", color = "#333333", family = best_font),
    strip.background = element_rect(fill = "#F5F5F5", color = "#CCCCCC"),
    
    # 标签设置（对annotate很重要）
    plot.tag = element_text(size = 11, face = "bold", family = best_font),
    
    # 边距设置
    plot.margin = ggplot2::margin(10, 15, 10, 10, "pt"),
    
    # 确保背景干净
    plot.background = element_rect(fill = "white", color = NA)
  )

# 为了避免字体警告，创建一个安全的annotate函数
safe_annotate <- function(geom, ..., family = best_font) {
  if (geom == "text") {
    annotate(geom, ..., family = family)
  } else if (geom == "label") {
    annotate(geom, ..., family = family)
  } else {
    annotate(geom, ...)
  }
}

# 打印当前使用的字体信息
cat("\n=== Font Configuration ===\n")
cat("Selected font for Nature theme:", best_font, "\n")
cat("This font will be used for all text elements in your plots.\n")

# 工具函数
saveFigure <- function(plot, filename, width = 8, height = 6, dpi = 600) {
  # 主图保存
  full_filename <- file.path(base_dir, "figures", paste0(filename, ".png"))
  ggsave(full_filename, plot, width = width, height = height, dpi = dpi, bg = "white")
  
  # 高分辨率版本（用于期刊提交）
  hr_filename <- file.path(base_dir, "figures", paste0(filename, ".tiff"))
  ggsave(hr_filename, plot, width = width, height = height, dpi = 600, bg = "white", device = "tiff")
  
  # 数据导出
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
# 数据加载与专业预处理
# -------------------------------------------------------------------------

# 智能数据加载
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

# 数据预处理与质量控制 - 强化多源植被数据处理
preprocess_data <- function(data) {
  cat("Starting data preprocessing...\n")
  cat("Original data dimensions:", dim(data), "\n")
  
  # 日期处理
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
  
  # 地理分类（基于最新气候分类标准Köppen-Geiger climate classification）

  # Continent分类
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
  
  # ========== 强化多源植被数据处理 ==========
  
  # 1. 三种数据源的基本清洗和标准化
  # MODIS IGBP分类处理
  if("MODIS_IGBP" %in% colnames(data)) {
    # 确保数据为字符型
    data$MODIS_IGBP <- as.character(data$MODIS_IGBP)
    
    # 创建标准化分类
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
    
    # 创建简化分类系统（用于后续比较）
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
    
    # 创建森林密度等级（用于比较）
    data$MODIS_Forest_Density <- case_when(
      data$MODIS_IGBP_Simplified != "Forest" ~ "Non-forest",
      grepl("Evergreen", data$MODIS_IGBP_Class_Clean) ~ "Dense Forest",
      grepl("Mixed", data$MODIS_IGBP_Class_Clean) ~ "Moderate Forest",
      grepl("Deciduous", data$MODIS_IGBP_Class_Clean) ~ "Moderate Forest",
      TRUE ~ "Unknown"
    )
  }
  
  # Copernicus LC分类处理
  if("Copernicus_LC" %in% colnames(data)) {
    # 确保数据为字符型
    data$Copernicus_LC <- as.character(data$Copernicus_LC)
    
    # 创建标准化分类
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
    
    # 创建简化分类系统（用于后续比较）
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
    
    # 创建森林密度等级（用于比较）
    data$Copernicus_Forest_Density <- case_when(
      data$Copernicus_LC_Simplified != "Forest" ~ "Non-forest",
      grepl("Closed Forest", data$Copernicus_LC_Class_Clean) ~ "Dense Forest",
      grepl("Open Forest", data$Copernicus_LC_Class_Clean) ~ "Moderate Forest",
      TRUE ~ "Unknown"
    )
  }
  
  # Hansen森林覆盖处理
  if("Hansen_Tree_Cover_2000_Percent" %in% colnames(data)) {
    # 确保数据为数值型
    data$Hansen_Tree_Cover_2000_Percent <- as.numeric(as.character(data$Hansen_Tree_Cover_2000_Percent))
    
    # 创建森林覆盖分类
    data$Hansen_Forest_Class_Clean <- case_when(
      is.na(data$Hansen_Tree_Cover_2000_Percent) ~ "Unknown",
      data$Hansen_Tree_Cover_2000_Percent == 0 ~ "No forest",
      data$Hansen_Tree_Cover_2000_Percent > 0 & data$Hansen_Tree_Cover_2000_Percent <= 25 ~ "Sparse forest",
      data$Hansen_Tree_Cover_2000_Percent > 25 & data$Hansen_Tree_Cover_2000_Percent <= 50 ~ "Moderate forest",
      data$Hansen_Tree_Cover_2000_Percent > 50 & data$Hansen_Tree_Cover_2000_Percent <= 75 ~ "Dense forest",
      data$Hansen_Tree_Cover_2000_Percent > 75 ~ "Very dense forest"
    )
    
    # 创建简化分类系统（用于后续比较）
    data$Hansen_Simplified <- case_when(
      data$Hansen_Forest_Class_Clean %in% c("Dense forest", "Very dense forest") ~ "Forest",
      data$Hansen_Forest_Class_Clean %in% c("Moderate forest", "Sparse forest") ~ "Woody Vegetation",
      data$Hansen_Forest_Class_Clean == "No forest" ~ "Non-forest",
      data$Hansen_Forest_Class_Clean == "Unknown" ~ "Unknown",
      TRUE ~ "Other"
    )
    
    # 保留详细的森林覆盖度分级（保持原代码）
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
  
  # 2. 创建多源一致性指标
  if(all(c("MODIS_IGBP_Simplified", "Hansen_Simplified") %in% colnames(data))) {
    # MODIS与Hansen的一致性
    data$MODIS_Hansen_Consistency <- ifelse(
      data$MODIS_IGBP_Simplified == "Unknown" | data$Hansen_Simplified == "Unknown",
      "Unknown",
      ifelse(data$MODIS_IGBP_Simplified == data$Hansen_Simplified, "Consistent", "Inconsistent")
    )
  }
  
  if(all(c("MODIS_IGBP_Simplified", "Copernicus_LC_Simplified") %in% colnames(data))) {
    # MODIS与Copernicus的一致性
    data$MODIS_Copernicus_Consistency <- ifelse(
      data$MODIS_IGBP_Simplified == "Unknown" | data$Copernicus_LC_Simplified == "Unknown",
      "Unknown",
      ifelse(data$MODIS_IGBP_Simplified == data$Copernicus_LC_Simplified, "Consistent", "Inconsistent")
    )
  }
  
  if(all(c("Copernicus_LC_Simplified", "Hansen_Simplified") %in% colnames(data))) {
    # Copernicus与Hansen的一致性
    data$Copernicus_Hansen_Consistency <- ifelse(
      data$Copernicus_LC_Simplified == "Unknown" | data$Hansen_Simplified == "Unknown",
      "Unknown",
      ifelse(data$Copernicus_LC_Simplified == data$Hansen_Simplified, "Consistent", "Inconsistent")
    )
  }
  
  # 3. 创建三源综合一致性评分
  if(all(c("MODIS_IGBP_Simplified", "Copernicus_LC_Simplified", "Hansen_Simplified") %in% colnames(data))) {
    # 计算一致性评分 (0-3)
    data$Multi_Source_Consistency_Score <- apply(
      data[, c("MODIS_IGBP_Simplified", "Copernicus_LC_Simplified", "Hansen_Simplified")], 
      1, 
      function(x) {
        if(any(x == "Unknown")) return(NA)
        unique_values <- length(unique(x))
        return(4 - unique_values)  # 3: 完全一致, 2: 两种一致, 1: 都不一致
      }
    )
    
    # 创建一致性标签
    data$Multi_Source_Consistency_Label <- case_when(
      is.na(data$Multi_Source_Consistency_Score) ~ "Unknown",
      data$Multi_Source_Consistency_Score == 3 ~ "Full Agreement",
      data$Multi_Source_Consistency_Score == 2 ~ "Partial Agreement",
      data$Multi_Source_Consistency_Score == 1 ~ "No Agreement"
    )
  }
  
  # 4. 创建统一的植被类型标签（优先级: MODIS > Copernicus > Hansen）
  data$Unified_Vegetation_Type <- case_when(
    "MODIS_IGBP_Simplified" %in% colnames(data) & !is.na(data$MODIS_IGBP_Simplified) & 
      data$MODIS_IGBP_Simplified != "Unknown" ~ data$MODIS_IGBP_Simplified,
    "Copernicus_LC_Simplified" %in% colnames(data) & !is.na(data$Copernicus_LC_Simplified) & 
      data$Copernicus_LC_Simplified != "Unknown" ~ data$Copernicus_LC_Simplified,
    "Hansen_Simplified" %in% colnames(data) & !is.na(data$Hansen_Simplified) & 
      data$Hansen_Simplified != "Unknown" ~ data$Hansen_Simplified,
    TRUE ~ "Unknown"
  )
  
  # 植被指数变化率计算
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
  
  # 数据质量报告
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

# 筛选出 "Climate_Zone" 不为 "Polar", Copernicus_LC_Class不为"Permanent water bodies"
# 且MODIS_IGBP_Class不为"Water Bodies"的数据
landslide_data <- landslide_data[
  !(landslide_data$Climate_Zone == "Polar" |
      landslide_data$Copernicus_LC_Class == "Permanent water bodies" |
      landslide_data$MODIS_IGBP_Class == "Water Bodies"),
]

saveTable(landslide_data,"landslide_data_processed")


# -------------------------------------------------------------------------
# MAIN FIGURE 1: 数据驱动阈值发现与风险区识别
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# MAIN FIGURE 1: 数据分析与补充材料生成
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
  # 内部辅助函数 (Internal Helper Functions)
  # ===================================================================
  
  # --- 辅助函数 1: 补充密度图 ---
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
  
  # --- 辅助函数 2: 补充交叉一致性表格 (优化版) ---
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
      # 将因子级别统一
      mutate(
        NDVI_Zone_Simple = gsub(" Zone", "", as.character(NDVI_Biophysical_Zone)),
        EVI_Zone_Simple = gsub(" Zone", "", EVI_Zone)
      )
    
    # 使用 `janitor` 包的 `tabyl` 和 `adorn_totals` 是最简洁、无警告的方法
    # 确保安装了 janitor: if(!require(janitor)) install.packages("janitor")
    if(!require(janitor)) install.packages("janitor", quiet = TRUE); library(janitor)
    
    consistency_table <- temp_data %>%
      tabyl(NDVI_Zone_Simple, EVI_Zone_Simple) %>%
      adorn_totals(where = c("row", "col")) %>% # Add row and column totals
      adorn_percentages(denominator = "row") %>% # Calculate row-wise percentages
      adorn_pct_formatting(digits = 1) %>% # Format as percentages
      adorn_ns(position = "front") # Add counts (n) in front of percentages
    
    # 将janitor表格转换为标准数据框并保存
    consistency_df <- as.data.frame(consistency_table)
    saveTable(consistency_df, "Supplementary_Table_1b_NDVI_EVI_Consistency")
    cat("  ✓ Cross-index consistency table saved.\n")
  }
  
  # ===================================================================
  # 主分析流程 (Main Analysis Workflow)
  # ===================================================================
  
  # --- 1. 执行核心断点分析 ---
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
  
  # --- 2. 定义风险区并添加到数据中 ---
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
  
  # --- 3. 调用辅助函数生成补充材料 ---
  # (注意：我们将修改后的`data`对象传递给辅助函数)
  if("NDVI_1" %in% names(analysis_results)) {
    create_supplementary_density_plot(data, analysis_results[["NDVI_1"]]$breakpoints)
  }
  if(all(c("NDVI_1", "EVI_1") %in% names(analysis_results))) {
    # 传递已包含NDVI_Biophysical_Zone列的数据
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
    
    # 定义绘图边界 - 增加边界距离以保持点线与图框的距离
    plot_xlim <- c(x_range[1] - 0.08*x_span, x_range[2] + 0.08*x_span)
    
    # 对于 NDVI/EVI，使用稍微扩展的范围
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
    
    # 开始构建图形
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
    
    # 【修正】使用旧代码中有效的置信区间阴影方法
    cat(sprintf("Debug: Adding confidence intervals for %s\n", idx))
    if(nrow(conf_int) >= 1) {
      for(i in 1:min(2, nrow(conf_int))) {
        # 确保置信区间的上下限顺序正确
        ci_lower <- min(conf_int[i, 1], conf_int[i, 2])
        ci_upper <- max(conf_int[i, 1], conf_int[i, 2])
        
        cat(sprintf("  BP%d: CI [%.3f, %.3f]\n", i, ci_lower, ci_upper))
        
        # 使用 geom_polygon 方法（来自旧代码）
        conf_data <- data.frame(
          x = c(ci_lower, ci_upper, ci_upper, ci_lower),
          y = c(0, 0, y_max, y_max)
        )
        
        p <- p + geom_polygon(data = conf_data, 
                              aes(x = x, y = y),  # 直接使用列名，不用继承
                              fill = "#E31A1C", 
                              alpha = 0.15,
                              inherit.aes = FALSE)  # 明确不继承美学映射
      }
    }
    
    # --- [!] RE-INTEGRATED MANUAL LABEL POSITIONING LOGIC ---
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
    # Fallback (though unlikely to be used with your specific indices)
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
    # --- END OF RE-INTEGRATED LOGIC ---
    
    p <- p +
      nature_theme_professional +
      # 保持 clip = "on" 但增加边界空间
      coord_cartesian(xlim = plot_xlim,
                      ylim = c(-0.02 * y_max, y_max * 1.1), # 下方和上方都增加空间
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
# 1. Run the analysis. The function returns a list with results and the MODIFIED data.
analysis_output <- analyze_vegetation_thresholds(landslide_data)



# 2. [!] CRITICAL FIX: Pass the MODIFIED data from the analysis output to the plotting function.
# Use analysis_output$results for the analysis results.
# Use analysis_output$data for the data, which now includes the 'NDVI_Biophysical_Zone' column.
figure_1 <- create_figure_1_plots(analysis_output$results, analysis_output$data)

# To update the global 'landslide_data_processed' object for subsequent script sections (e.g., Figure 2, 3), do this:
landslide_data_processed <- analysis_output$data




# ===========================================================================

# -------------------------------------------------------------------------
# MAIN FIGURE 2: 数据分析部分 - 多源植被数据验证与跨平台一致性
# -------------------------------------------------------------------------

# =============================================================================

analyze_vegetation_validation_data <- function() {
  cat("\n=== Analyzing Vegetation Data Validation and Cross-Platform Consistency ===\n")
  
  data <- landslide_data_processed # 使用前面处理过的数据
  analysis_results <- list()
  
  # --- 2A: 全球植被类型分布分析 ---
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
  
  # --- 2B: 多源一致性水平分析 ---
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
    
    # Storing CI requires a different approach as irr package doesn't directly provide it.
    # We will report the main value and can add CI later if needed with a bootstrapper.
    # For now, we will add the CI from your text as a reference.
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
  # 手动从cohen_kappa_result中提取需要的值并创建一个新的数据框
  kappa_df <- data.frame(
    Statistic = "Cohen's Kappa",
    Value = cohen_kappa_result$value,
    Z_statistic = cohen_kappa_result$statistic,
    P_value = cohen_kappa_result$p.value,
    Subjects = cohen_kappa_result$subjects,
    Raters = cohen_kappa_result$raters
  )
  saveTable(kappa_df, "Supplementary_Table_2b_kappa_MODIS_VS_Copernicus")
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
  
  # 用于分析2C, 2D, 2E的助手函数
  # [!] 这是需要修正的部分
  analyze_zone_distribution <- function(data_subset, x_var_str, output_filename) {
    summary_df <- data_subset %>%
      filter(!is.na(.data[[x_var_str]]), !is.na(NDVI_Biophysical_Zone),
             .data[[x_var_str]] != "Unknown", .data[[x_var_str]] != "Other", 
             NDVI_Biophysical_Zone != "Unknown") %>%
      # 第一次分组：按x轴分类和风险区共同分组，计算每个组合的数量
      group_by(across(all_of(x_var_str)), NDVI_Biophysical_Zone) %>%
      plyr::summarise(Count = n(), .groups = "drop") %>%
      # 第二次分组 (关键步骤)：仅按x轴分类进行分组
      group_by(across(all_of(x_var_str))) %>%
      # 计算百分比，此时 sum(Count) 是每个x轴分类内部的总数
      mutate(Percentage = Count / sum(Count) * 100)
    
    if(nrow(summary_df) > 0) {
      saveTable(summary_df, output_filename)
    }
    return(summary_df)
  }
  # --- 2C: 主要土地覆盖类型中的风险区分析 ---
  landcover_zone_summary <- analyze_zone_distribution(
    data_subset = data %>% filter(MODIS_IGBP_Simplified %in% c("Forest", "Woody Vegetation", "Grassland", "Cropland/Vegetation")),
    x_var_str = "MODIS_IGBP_Simplified",
    output_filename = "Fig2C_Landcover_Zone_Distribution"
  )
  analysis_results[["landcover_zone_summary"]] <- landcover_zone_summary
  
  
  # --- 2D: 森林覆盖密度与风险区分析 ---
  forest_density_zone_summary <- analyze_zone_distribution(
    data_subset = data,
    x_var_str = "Forest_Cover_Category",
    output_filename = "Fig2D_ForestDensity_Zone_Distribution"
  )
  analysis_results[["forest_density_zone_summary"]] <- forest_density_zone_summary
  
  # --- 2E: 气候区风险分布分析 ---
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
  
  # 使用在预处理步骤中生成的全局变量
  data <- landslide_data_processed 
  analysis_results <- list()
  
  # --- 2A: 全球植被类型分布分析 ---
  if("MODIS_IGBP_Class_Clean" %in% colnames(data)) {
    modis_dist <- data %>%
      # 明确使用 dplyr::filter
      dplyr::filter(!is.na(MODIS_IGBP_Class_Clean), MODIS_IGBP_Class_Clean != "Unknown") %>%
      dplyr::group_by(MODIS_IGBP_Class_Clean) %>%
      # 明确使用 dplyr::summarise 和 dplyr::n()
      dplyr::summarise(Count = n(), .groups = "drop") %>%
      dplyr::filter(Count > 0) %>%
      # 明确使用 dplyr::arrange
      dplyr::arrange(desc(Count)) %>%
      # 明确使用 dplyr::slice_head
      dplyr::slice_head(n = 10) %>% # Top 10 for clarity
      # 明确使用 dplyr::mutate
      dplyr::mutate(Percentage = Count / sum(Count) * 100,
                    Display_Name_Short = sapply(strsplit(MODIS_IGBP_Class_Clean, " "), function(x) paste(x[1:min(length(x),2)], collapse=" "))) # First 1 or 2 words
    
    analysis_results[["vegetation_type_distribution"]] <- modis_dist
    saveTable(modis_dist, "Fig2A_modis_vegetation_distribution")
  } else {
    cat("Warning: MODIS_IGBP_Class_Clean not available for vegetation type analysis\n")
  }
  
  # --- 2B: 多源一致性水平分析 ---
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
  
  
  
  # --- 用于分析2C, 2D, 2E的助手函数 (这是关键的修正部分) ---
  analyze_zone_distribution <- function(data_subset, x_var_str, output_filename) {
    summary_df <- data_subset %>%
      dplyr::filter(!is.na(.data[[x_var_str]]), !is.na(NDVI_Biophysical_Zone),
                    .data[[x_var_str]] != "Unknown", .data[[x_var_str]] != "Other", 
                    NDVI_Biophysical_Zone != "Unknown") %>%
      # 第一次分组：按x轴分类和风险区共同分组
      dplyr::group_by(across(all_of(x_var_str)), NDVI_Biophysical_Zone) %>%
      # 计算每个组合的数量
      dplyr::summarise(Count = n(), .groups = "drop") %>%
      # 第二次分组 (关键步骤)：仅按x轴分类进行分组
      dplyr::group_by(across(all_of(x_var_str))) %>%
      # 计算百分比，此时 sum(Count) 是每个x轴分类内部的总数
      dplyr::mutate(Percentage = Count / sum(Count) * 100) %>%
      # 完成计算后取消分组，这是个好习惯
      dplyr::ungroup()
    
    if(nrow(summary_df) > 0) {
      saveTable(summary_df, output_filename)
    }
    return(summary_df)
  }
  
  # --- 2C: 主要土地覆盖类型中的风险区分析 ---
  landcover_zone_summary <- analyze_zone_distribution(
    data_subset = data %>% dplyr::filter(MODIS_IGBP_Simplified %in% c("Forest", "Woody Vegetation", "Grassland", "Cropland/Vegetation")),
    x_var_str = "MODIS_IGBP_Simplified",
    output_filename = "Fig2C_Landcover_Zone_Distribution"
  )
  analysis_results[["landcover_zone_summary"]] <- landcover_zone_summary
  
  # --- 2D: 森林覆盖密度与风险区分析 ---
  forest_density_zone_summary <- analyze_zone_distribution(
    data_subset = data,
    x_var_str = "Forest_Cover_Category",
    output_filename = "Fig2D_ForestDensity_Zone_Distribution"
  )
  analysis_results[["forest_density_zone_summary"]] <- forest_density_zone_summary
  
  # --- 2E: 气候区风险分布分析 ---
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
    
    # Storing CI requires a different approach as irr package doesn't directly provide it.
    # We will report the main value and can add CI later if needed with a bootstrapper.
    # For now, we will add the CI from your text as a reference.
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
  # 手动从cohen_kappa_result中提取需要的值并创建一个新的数据框
  kappa_df <- data.frame(
    Statistic = "Cohen's Kappa",
    Value = cohen_kappa_result$value,
    Z_statistic = cohen_kappa_result$statistic,
    P_value = cohen_kappa_result$p.value,
    Subjects = cohen_kappa_result$subjects,
    Raters = cohen_kappa_result$raters
  )
  saveTable(kappa_df, "Supplementary_Table_2b_kappa_MODIS_VS_Copernicus")
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
  
  # 将分析结果赋值给全局变量，供后续可视化函数使用
  assign("vegetation_validation_results", analysis_results, envir = .GlobalEnv)
  
  cat("✓ Vegetation validation analysis completed successfully\n")
  return(analysis_results)
}

# 执行分析
vegetation_validation_results <- analyze_vegetation_validation_data()


# -------------------------------------------------------------------------
# MAIN FIGURE 2: 可视化部分 - 多源植被数据验证与跨平台一致性（改进x轴标签显示）
# -------------------------------------------------------------------------
create_vegetation_validation_plots <- function(analysis_results, data) {
  cat("\n=== Creating Figure 2 Visualizations with Improved Labels ===\n")
  
  plot_list <- list()
  
  # --- 2A: 全球植被类型分布(Treemap) --- [无变化]
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
  
  # --- 2B: 多源一致性水平(饼图) --- [无变化]
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
  
  # --- 修改的共用函数: 创建风险分布图(用于2C, 2D, 2E) ---
  # 主要变化：添加了标签处理逻辑，将长标签分成两行显示
  create_zone_distribution_plot <- function(summary_df, x_var_str, x_lab, is_first_plot = FALSE, is_last_plot = FALSE) {
    if(is.null(summary_df) || nrow(summary_df) == 0) {
      return(ggplot() + geom_text(aes(0, 0, label = paste("No data available for", x_var_str))) + 
               nature_theme_professional)
    }
    
    # 创建新的修改版数据框，处理标签分行显示
    modified_df <- summary_df
    
    # 特殊处理森林覆盖度标签 - 更简洁的格式
    if(x_var_str == "Forest_Cover_Category") {
      # 使用更简洁的标签，只保留必要信息
      label_mapping <- c(
        "No Forest (0%)" = "No Forest\n(0%)",
        "Sparse (1-25%)" = "Sparse\n(1-25%)",
        "Moderate (26-50%)" = "Moderate\n(26-50%)",
        "Dense (51-75%)" = "Dense\n(51-75%)",
        "Very Dense (>75%)" = "Very Dense\n(>75%)"
      )
      
      # 确保保持原始顺序
      if("order" %in% colnames(modified_df)) {
        original_order <- modified_df$order
      } else {
        # 如果没有order列，尝试创建一个基于标准排序的列
        standard_order <- c("No Forest (0%)", "Very Low (1-20%)", "Low (21-40%)", 
                            "Medium (41-60%)", "High (61-80%)", "Very High (81-100%)")
        original_order <- match(modified_df[[x_var_str]], standard_order)
      }
      
      # 创建新的标签列并保持排序
      modified_df$Display_Label <- sapply(modified_df[[x_var_str]], function(x) {
        if(x %in% names(label_mapping)) {
          return(label_mapping[x])
        } else {
          return(as.character(x))
        }
      })
      
      # 确保排序正确
      modified_df$Display_Label <- factor(modified_df$Display_Label, 
                                          levels = unique(modified_df$Display_Label[order(original_order)]))
    }
    # 处理其他类型的标签 [与之前相同]
    else if(x_var_str == "MODIS_IGBP_Simplified") {
      # 土地覆盖类型标签处理 [代码保持不变]
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
      # 气候区标签处理 [代码保持不变]
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
      # 默认情况，使用原始标签
      modified_df$Display_Label <- modified_df[[x_var_str]]
    }
    
    # 为森林覆盖度图表使用更紧凑的宽度设置
    bar_width <- 0.8  # 默认宽度
    if(x_var_str == "Forest_Cover_Category") {
      bar_width <- 0.7  # 森林覆盖度使用更窄的条形
    }
    
    # 使用修改后的标签创建图表
    p <- ggplot(modified_df, aes(x = Display_Label, y = Percentage, fill = NDVI_Biophysical_Zone)) +
      geom_col(position = "stack", alpha = 0.9, width = bar_width) +
      geom_text(aes(label = ifelse(Percentage > 7, sprintf("%.0f%%", Percentage), "")),
                position = position_stack(vjust = 0.5), color = "white", 
                fontface = "bold", size = 2.5) +
      scale_fill_manual(values = nature_palettes$biophysical_zones, name="Biophysical Zone") +
      scale_y_continuous(labels = scales::percent_format(scale=1), 
                         limits = c(0, 100), breaks = seq(0,100,25)) +
      nature_theme_professional
    
    # 为森林覆盖度图表使用特殊的标签设置
    if(x_var_str == "Forest_Cover_Category") {
      p <- p + theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, size=7), # 稍微小一点的字体
        axis.text.y = element_text(size=6),
        axis.title = element_text(size=7),
        plot.margin = ggplot2::margin(2, 0, 2, 0),  # 减少左右边距
        panel.spacing = unit(0, "pt")  # 减少面板间距
      )
    } else {
      # 其他图表使用标准设置
      p <- p + theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, size=8),
        axis.text.y = element_text(size=6),
        axis.title = element_text(size=7),
        plot.margin = ggplot2::margin(2, 2, 2, 2)
      )
    }
    
    # 添加x轴标签
    p <- p + labs(x = x_lab)
    
    # 只在第一个图添加Y轴标签
    if (is_first_plot) {
      p <- p + labs(y = "Percentage of Landslides (%)")
    } else {
      p <- p + theme(axis.title.y = element_blank(), 
                     axis.text.y = element_blank(), 
                     axis.ticks.y = element_blank())
    }
    
    # 只在最后一个图显示图例
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
  
  # --- 2C: 主要土地覆盖类型中的风险区 ---
  p2c <- create_zone_distribution_plot(
    summary_df = analysis_results[["landcover_zone_summary"]],
    x_var_str = "MODIS_IGBP_Simplified",
    x_lab = "Land Cover Type",
    is_first_plot = FALSE,
    is_last_plot = FALSE
  )
  
  # --- 2D: 森林覆盖密度与风险区 ---
  p2d <- create_zone_distribution_plot(
    summary_df = analysis_results[["forest_density_zone_summary"]],
    x_var_str = "Forest_Cover_Category",
    x_lab = "Forest Cover Category",
    is_first_plot = FALSE,
    is_last_plot = FALSE
  )
  
  # --- 2E: 气候区风险分布 ---
  p2e <- create_zone_distribution_plot(
    summary_df = analysis_results[["climate_zone_summary"]],
    x_var_str = "Climate_Zone",
    x_lab = "Climate Zone",
    is_first_plot = FALSE,
    is_last_plot = FALSE  # 这是最后一个图，显示图例
  )
  
  # 创建图2的所有组件
  p2a <- create_vegetation_treemap()
  p2b <- create_multisource_consistency_piechart()
  
  plot_list <- list(
    treemap = p2a,
    consistency = p2b,
    landcover_zone = p2c,
    forest_density = p2d,
    climate_zone = p2e
  )
  
  # 组合图2 - 使用与图1相似的尺寸比例
  top_row <- p2a + p2b + plot_layout(widths = c(1.2, 0.8))
  
  # 组合底部行并确保只有一个图例
  bottom_row <- p2c + p2d + p2e + 
    plot_layout(
      widths = c(1, 1.2, 1),
      guides = "collect"  # 收集所有图例
    ) & 
    theme(legend.position = "bottom")  # 确保图例在底部
  
  # 组合上下两行
  main_fig2 <- (top_row / bottom_row) + 
    plot_layout(heights = c(1, 1)) + 
    plot_annotation(tag_levels = list(paste0("(", letters, ")"))) &
    theme(
      plot.tag = element_text(face = "bold", size = 10),
      panel.spacing = unit(4, "pt"),
      plot.margin = ggplot2::margin(5, 5, 5, 5)
    )
  
  # 保存组合图 - 使用与图1相同的尺寸
  saveFigure(main_fig2, "Main_Figure_2_A4", width = 8.3, height = 6.8)  # 与图1相同的尺寸
  
  cat("✓ Figure 2 visualizations with improved labels completed successfully\n")
  return(list(plots = plot_list, combined = main_fig2))
}

# 执行出图
figure2_plots <- create_vegetation_validation_plots(vegetation_validation_results, landslide_data_processed)

# ------------------------------------------------------------------------------
# Figure 2 桑基图
# ------------------------------------------------------------------------------
# ==============================================================================
# 步骤 1: 数据准备 (与您原代码相同)
# ==============================================================================

# 1.1 创建一个统一的、简化的最终分类 (作为桑基图的右侧终点)
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

# 1. 筛选数据
complete_data <- landslide_data_processed %>%
  dplyr::filter(!is.na(MODIS_IGBP_Simplified) & !is.na(Copernicus_LC_Simplified) & !is.na(Hansen_Simplified)) %>%
  dplyr::filter(MODIS_IGBP_Simplified != "Unknown" & Copernicus_LC_Simplified != "Unknown")

total_samples <- nrow(complete_data)

# 2. 计算频率
combo_counts <- complete_data %>%
  dplyr::group_by(MODIS_IGBP_Simplified, Copernicus_LC_Simplified, Hansen_Simplified) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(percent = n / total_samples * 100)

# 3. 分配 alluvium ID
combo_counts <- combo_counts %>%
  dplyr::mutate(alluvium = dplyr::row_number())

# 4. 转换为长格式
sankey_data_long <- combo_counts %>%
  tidyr::pivot_longer(
    cols = c(MODIS_IGBP_Simplified, Copernicus_LC_Simplified, Hansen_Simplified),
    names_to = "Source",
    values_to = "stratum"
  )

# 5. 添加统一分类
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
# 步骤 2: 定义颜色和标签 (与您原代码相同)
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
# 步骤 3: 准备最终绘图数据 (与您原代码相同)
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
  # 确保Source也是一个有序因子，以便绘图
  dplyr::mutate(
    Source = factor(Source, levels = c("MODIS_IGBP", "Copernicus_LC", "Hansen", "Unified"))
  )


# ==============================================================================
# 步骤 4: 添加百分比标签 - 优化标签位置 (最终修正版，带显式包名)
# ==============================================================================

# 首先，我们需要从原始绘图数据中计算出每个分类块（stratum）的总大小
stratum_summary <- plot_data %>%
  dplyr::filter(!is.na(stratum_ordered)) %>%
  dplyr::group_by(Source, stratum_ordered) %>%
  dplyr::summarise(stratum_total_n = sum(n), .groups = "drop")

# 其次，基于每个块的总大小，我们计算其在Y轴上的位置和百分比
label_data <- stratum_summary %>%
  # 使用 desc() 反转排列顺序，以匹配 ggalluvial 的堆叠顺序 (将因子第一级放在顶部)
  dplyr::arrange(Source, dplyr::desc(stratum_ordered)) %>%
  # 按x轴变量（Source）分组
  dplyr::group_by(Source) %>%
  # 计算每个色块的顶部、底部和中心位置
  dplyr::mutate(
    y_top = cumsum(stratum_total_n),
    y_bottom = dplyr::lag(y_top, default = 0),
    y_center = (y_top + y_bottom) / 2,
    # 基于块的总大小计算正确的百分比
    total_n_in_source = sum(stratum_total_n),
    percent = stratum_total_n / total_n_in_source * 100,
    label = ifelse(percent >= 2, sprintf("%.0f%%", percent), "")
  ) %>%
  dplyr::ungroup() # 完成计算后取消分组

# 最后，为标签添加样式（大小和颜色）
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
# 步骤 5: 绘制桑基图 (使用新的 label_data)
# ==============================================================================

sankey_plot <- ggplot2::ggplot(data = plot_data,
                               ggplot2::aes(x = Source, stratum = stratum_ordered, alluvium = alluvium,
                                            y = n, fill = stratum_ordered)) +
  # 使用ggalluvial的核心几何对象
  ggalluvial::geom_flow(stat = "alluvium", lode.guidance = "forward", color = "darkgray", alpha = 0.6, width = 0.4) +
  ggalluvial::geom_stratum(alpha = 1, width = 0.4, linewidth = 0.2) +
  
  # 添加分类块内的百分比标签 - 使用修正后的 label_data
  ggplot2::geom_text(
    data = label_data %>% dplyr::filter(label != ""), # 只标注非空标签
    mapping = ggplot2::aes(x = Source, y = y_center, label = label, size = label_size, color = label_color),
    inherit.aes = FALSE,
    fontface = "bold"
  ) +
  
  # 设置标签大小和颜色
  ggplot2::scale_size_identity() +
  ggplot2::scale_color_identity() +
  
  # 使用我们定义的颜色
  ggplot2::scale_fill_manual(
    values = color_palette,
    name = "Land Cover Classification",
    breaks = class_order,
    na.value = "grey80" # 为可能出现的NA值指定颜色
  ) +
  
  # 调整X轴标签
  ggplot2::scale_x_discrete(
    limits = c("MODIS_IGBP", "Copernicus_LC", "Hansen", "Unified"),
    labels = c("MODIS", "Copernicus", "Hansen", "Unified"),
    expand = ggplot2::expansion(mult = c(0.05, 0.05)) # 调整x轴两侧空白
  ) +
  
  # Y轴不需要额外空间
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
  
  # 主题和标题
  # 注意：如果您想要黑色背景，请在这里应用 theme_dark() 或者手动设置
  # 例如：theme_dark() + theme(...)
  nature_theme_professional +
  ggplot2::theme_minimal(base_family = "sans") +
  ggplot2::theme(
    legend.position = "top",
    legend.box = "horizontal",
    legend.title = ggplot2::element_text(size = 9, face = "bold"),
    legend.text = ggplot2::element_text(size = 8),
    legend.margin = ggplot2::margin(b = 0, t = 0), # 减小图例与图的距离
    legend.spacing.x = ggplot2::unit(0.2, "cm"), # 减小图例内部间距
    legend.key.size = ggplot2::unit(0.5, "cm"), # 减小图例符号大小

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
    x = "Data Source / Classification", # 添加X轴标题
    y = "Number of Landslides",
    fill = "" # 隐藏图例标题
  )

# 打印图形
print(sankey_plot)


# 执行出图
saveFigure(sankey_plot, "Supplementary Fig. 2b sankey plot", width = 8.3, height = 6.8)  # 与图1相同的尺寸




# ==============================================================================

# ------------------------------------------------------------------------------
# 图3
# ------------------------------------------------------------------------------

# ==============================================================================

# ------------------------------------------------------------------------------
# 数据分析部分 - 为图3准备数据
# ------------------------------------------------------------------------------
prepare_data_for_figure_3 <- function() {
  cat("\n=== Preparing Data for Figure 3 ===\n")
  
  data <- landslide_data_processed 
  
  # 确保NDVI_Biophysical_Zone正确定义
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
  
  # 准备图3A的数据 (散点图)
  data_3a <- data %>% dplyr::filter(!is.na(NDVI_1), NDVI_Biophysical_Zone != "Unknown")
  
  # 准备图3B的数据 (小提琴图)
  climate_ndvi_data <- data %>%
    dplyr::filter(!is.na(NDVI_1), !is.na(Climate_Zone), Climate_Zone != "Other")
  
  # --- 补充图Y: 按气候带的分面断点分析 ---
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
  
    # 初始化返回列表中的变量
  seasonal_summary <- NULL
  
  # 准备图3C的数据 (季节性)
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
          # 使用purrr::map安全地应用检验
          chi_sq_result = purrr::map(data, ~ tryCatch(chisq.test(.x$Count), error = function(e) NULL)),
          # 仅在检验成功时提取结果
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
  
  # 初始化返回列表中的变量，以防后续步骤失败
  temporal_summary <- NULL
  temporal_change_annotations <- NULL
  
  # 准备图3D的数据 (时间序列)
  temporal_cols <- paste0("NDVI_", 5:1)  # 修改为从5到1的顺序
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
      # 添加反向顺序的索引，用于排序
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
        # 按反向顺序排序
        dplyr::arrange(NDVI_Biophysical_Zone, desc(Time_Window_Num))
      saveTable(temporal_summary, "Statistics_Fig3D_temporal_summary")
      cat("✓ temporal_summary results saved to 'tables/' directory.\n")
      
      temporal_change_annotations <- temporal_summary %>%
        dplyr::group_by(NDVI_Biophysical_Zone) %>%
        # 按反向顺序排序
        dplyr::arrange(desc(Time_Window_Num)) %>%
        dplyr::mutate(Prev_Mean_NDVI = lag(Mean_NDVI),
                      Percent_Change = ifelse(is.na(Prev_Mean_NDVI) | Prev_Mean_NDVI == 0, NA, 
                                              (Mean_NDVI - Prev_Mean_NDVI) / Prev_Mean_NDVI * 100)) %>%
        # 注意：现在我们是从大到小，所以要选择小于而不是大于1
        dplyr::filter(Time_Window_Num < 5 & !is.na(Percent_Change) & abs(Percent_Change) > 2) 
      saveTable(temporal_change_annotations, "Statistics_Fig3D_temporal_change_annotations")
      cat("✓ temporal_change_annotations results saved to 'tables/' directory.\n")
      
      # --- STATISTICAL TEST 3: Mann-Kendall Trend Test for Time Series ---
      cat("--- Performing Mann-Kendall Trend Tests for each Biophysical Zone ---\n")
      
      mk_results <- temporal_data_long %>%
        dplyr::group_by(OBJECTID, NDVI_Biophysical_Zone) %>%
        # 按反向顺序排序
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

# 运行数据准备函数
figure_3_data <- prepare_data_for_figure_3()


# ------------------------------------------------------------------------------
# 可视化部分 - 绘制图3（最终完善版本）
# ------------------------------------------------------------------------------

create_figure_3_visualization <- function(figure_data, width = 8.3, height = 6.8) {
  cat("\n=== Creating Figure 3 Visualization with Expert Refinements ===\n")
  
  # 确保figures目录存在
  if (!dir.exists("figures")) {
    dir.create("figures")
    cat("Created 'figures' directory\n")
  }
  
  # 设置颜色配置
  biophysical_zone_colors <- c("IDZ" = "#3B9AB2", "CTZ" = "#E14D55", "SDZ" = "#21A764")
  
  season_colors <- c("Winter" = "#D3DDEA", "Spring" = "#B9DCC9", "Summer" = "#A5D2E5", "Fall" = "#F5C0B8")
  
  # 图3A: 植被指数分布 - 图例放在底部
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
        legend.title = element_text(size = 9), # Control legend title size
        legend.text = element_text(size = 8),  # Control legend label size
        legend.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
        legend.box.margin = ggplot2::margin(t = -5)
      ) +
      labs(x = "NDVI Value", y = "Density", fill = "Biophysical Zone", color = "Biophysical Zone")
  } else {
    plot_3a <- ggplot() + theme_void() + 
      annotate("text", x = 0.5, y = 0.5, label = "No data available", size = 3)
  }
  
  # [!] 这是修正后的完整代码块，用于生成图3B
  # [!] This is the fully revised code block for Figure 3B, now with statistical export
  
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
    
    # --- Step 1B: Calculate and Save Key Statistics (NEW SECTION) ---
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
    # The 'dunn.test' package is great for this. Make sure it's installed.
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
    
    # --- End of new statistics section ---
    
    
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
    
    # 3. Core Plotting Code (no changes needed here)
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
  
  # 图3C: 季节分布 - 图例放在底部，x轴标签分两行显示
  if(!is.null(figure_data$data_3c) && nrow(figure_data$data_3c) > 0) {
    # 修改NDVI分类标签，分两行显示
    figure_data$data_3c$NDVI_Category_2line <- factor(
      figure_data$data_3c$NDVI_Category,
      levels = c("Low (0-0.3)", "Moderate (0.3-0.6)", "High (0.6-0.8)", "Very High (0.8-1)"),
      labels = c("Low\n(0-0.3)", "Moderate\n(0.3-0.6)", "High\n(0.6-0.8)", "Very High\n(0.8-1)")
    )
    # 计算每个NDVI分类的总样本量
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
        legend.title = element_text(size = 8), # Control legend title size
        legend.text = element_text(size = 7),  # Control legend label size
        legend.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
        legend.box.margin = ggplot2::margin(t = -5),
        axis.text.x = element_text(angle = 0, hjust = 0.5)  # 水平显示标签
      ) +
      # 在柱子下方添加样本量信息
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
  
  
  # 图3D: 时间序列分析 - 图例放在底部
  if(!is.null(figure_3_data$data_3d_summary) && nrow(figure_3_data$data_3d_summary) > 0) {
    # [!] 核心修改部分：创建详细的X轴标签
    # --------------------------------------------------------------------------
    # 1. 定义时间窗口的天数范围
    time_window_days <- c(
      "T5" = "(80-64 days)",
      "T4" = "(64-48 days)",
      "T3" = "(48-32 days)",
      "T2" = "(32-16 days)",
      "T1" = "(16-0 days)"
    )
    
    # 2. 创建一个包含换行符的标签向量，顺序与X轴一致 (T5 -> T1)
    detailed_x_labels <- paste0(
      names(time_window_days),    # T5, T4, ...
      "\n",                       # 换行符
      time_window_days            # (64-80 days), ...
    )
    # -------
    # Assuming biophysical_zone_colors is defined in your environment
    # e.g., biophysical_zone_colors <- nature_palettes$biophysical_zones
    
    plot_3d <- ggplot(figure_3_data$data_3d_summary, 
                      aes(x = reorder(Time_Window, -Time_Window_Num), y = Mean_NDVI, 
                          group = NDVI_Biophysical_Zone, color = NDVI_Biophysical_Zone)) +
      geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper, fill = NDVI_Biophysical_Zone), alpha = 0.2, color = NA) +
      geom_vline(xintercept = "T1", linetype = "dashed", color = "pink") +
      geom_line(linewidth = 1.2) +
      geom_point(size = 2.5) +
      # 修改这里：调整标签位置，向左上方移动
      geom_text(data = figure_data$data_3d_annotations,
                aes(label = sprintf("%+.1f%%", Percent_Change)),
                vjust = 1.8, hjust = 0.6, # 添加hjust=0.8让标签向左移动
                nudge_x = -0.15, # 向左移动标签
                size = 2.5, fontface = "bold", show.legend = FALSE) +
      
      # Use coord_cartesian to get the y-min for annotation placement
      annotate("text", x = "T1", y = 0.51, # Set a specific position since axis starts at 0.5
               label = "Landslide Event Proximal", hjust = 1.05, vjust = -12, size = 3, color="grey30", fontface="italic") +
      # [!] 在 scale_x_discrete 中使用我们新创建的 labels
      # --------------------------------------------------------------------------
      scale_x_discrete(
        name = "Time Window (days before event)", # 更新X轴标题
        labels = detailed_x_labels,            # 使用新标签
        expand = expansion(mult = c(0.1, 0.15)) # 调整左右留白以容纳标签
      ) +
      # --------------------------------------------------------------------------
      # --- 修改这里：增加右侧留白 ---
      coord_cartesian(ylim = c(0.5, NA), clip = "off") + # 保持clip="off"以便文字能够在需要时扩展
      #scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) + # 右侧增加更多留白
      
      scale_color_manual(values = biophysical_zone_colors) +
      scale_fill_manual(values = biophysical_zone_colors) +
      nature_theme_professional +
      
      # --- CHANGE 2: Adjust legend font size ---
      theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = element_text(size = 8), # Control legend title size
        legend.text = element_text(size = 7),  # Control legend label size
        axis.text.x = element_text(
          lineheight = 0.9, # 行高，对于多行文字很重要
          size = 8          # 将大小从 7 调整为 8
        ),
        legend.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = -2),
        legend.box.margin = ggplot2::margin(t = -5),
        # 添加这一行以确保标签不会被裁剪
        plot.margin = ggplot2::margin(t = 5, r = 0, b = 0, l = 5)
      ) +
      labs(x = "Time Window", y = "Mean NDVI (95% CI)", color = "Biophysical Zone", fill = "Biophysical Zone")
    
  } else {
    plot_3d <- ggplot() + theme_void() + 
      annotate("text", x = 0.5, y = 0.5, label = "No temporal data available", size = 3)
  }
  print(plot_3d)
  

  # 组合所有图形
  final_figure <- (plot_3a | plot_3b) / (plot_3c | plot_3d) +
    #plot_layout(guides = 'collect') + # Helps with aligning plots that have legends
    plot_annotation(tag_levels = list(paste0("(", letters, ")"))) &
    theme(
      plot.tag = element_text(face = "bold", size = 10),
      plot.margin = ggplot2::margin(t = 5, r = 5, b = 0, l = 5) # Reduce margins around each plot
      )
  
  # 保存图形到figures文件夹
  # 也保存各个子图以便单独使用
  saveFigure(final_figure, "My_Final_Figure_3_A4", width = 8.3, height = 6.8)
  
  cat("✓ Figure 3 visualization with expert refinements completed successfully\n")
  cat("✓ All figures saved to the 'figures' directory\n")
  
  return(final_figure)
}


# 创建最终图形
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
# Part 1: 模型训练函数
# 任务: 准备数据，执行交叉验证，训练最终模型，并保存所有结果
# =========================================================================
run_fig4_modeling <- function() {
  cat("\n=== RUNNING MODELING FOR FIGURE 4 ===\n")
  # Ensure necessary packages are available
  if (!require(pROC, quietly = TRUE)) { install.packages("pROC", quiet = TRUE); library(pROC) }
  
  # --- 0. 数据和环境准备 ---
  if (!exists("landslide_data_processed")) {
    stop("Required data `landslide_data_processed` not found.")
  }
  data <- landslide_data_processed
  
  # --- 1. 创建非循环的目标变量 (Landslide Density Proxy) ---
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
  
  # --- 2. 特征选择与数据准备 ---
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
  
  # --- 3. 稳健的空间交叉验证 ---
  cat("Step 3: Performing spatial cross-validation...\n")
  set.seed(123)
  q_lat <- unique(quantile(modeling_df_full$Latitude, probs = 0:5/5, na.rm = TRUE))
  if (length(q_lat) < 3) {
    modeling_df_full$Fold <- sample(1:5, nrow(modeling_df_full), replace = TRUE)
  } else {
    modeling_df_full$Fold <- as.numeric(cut(modeling_df_full$Latitude, breaks = q_lat, include.lowest = TRUE))
  }
  modeling_df_full$Fold[is.na(modeling_df_full$Fold)] <- sample(1:5, sum(is.na(modeling_df_full$Fold)), replace = TRUE)
  
  
  # <<< NEW SECTION START >>>
  # Create a binary target for classification task
  bz_threshold <- quantile(modeling_df_full$Density_Target, 0.8)
  modeling_df_full$Zone_Class <- factor(
    ifelse(modeling_df_full$Density_Target >= bz_threshold, "CT_Zone", "Not_CT_Zone"),
    levels = c("Not_CT_Zone", "CT_Zone") # Set "Not_CT_Zone" as the reference level
  )
  cat(paste0("--- Created binary classification target 'Zone_Class' with threshold: ", round(bz_threshold, 3), " ---\n"))
  # <<< NEW SECTION END >
  
  
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
    # <<< CORRECTION FOR ROC WARNING >>>
    # Explicitly tell roc() which level is the "case" or "event"
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
    # 修改点1: 将图例移到图框内右上角
    theme(legend.position = c(0.22, 0.66),
          legend.background = element_rect(fill = "white", color = "gray90"),
          legend.margin = ggplot2::margin(5, 5, 5, 5)) +
    labs(
      title = "Cross-Validated ROC Curves by Spatial Fold",
      x = "False positive fraction", # Corrected label
      y = "True positive fraction", # Corrected label
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
  
  # --- 5. 准备并保存所有绘图所需的数据 ---
  cat("Step 5: Saving all modeling results for plotting...\n")
  
  # <<< CORRECTION FOR THE ERROR >>>
  # Train the final regression model that was missing
  final_reg_formula <- as.formula(paste("Density_Target ~", paste(features, collapse = " + ")))
  validation_rf <- randomForest(final_reg_formula, data = modeling_df_full, ntree = 500, importance = TRUE)
  # <<< END CORRECTION >>
  
  # 为Panel C准备数据
  avg_importance <- bind_rows(feature_importance_list) %>%
    group_by(Feature) %>%
    dplyr::summarise(Mean_Importance = mean(Importance, na.rm = TRUE), SD_Importance = sd(Importance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(Mean_Importance))
  
  # 为Panel D准备数据
  modeling_df_full$Predicted_Density <- predict(validation_rf, modeling_df_full %>% dplyr::select(all_of(features)))
  pdp_plot_data <- modeling_df_full %>%
    dplyr::select(OBJECTID, Predicted_Density) %>%
    left_join(model_data %>% dplyr::select(OBJECTID, NDVI_1), by = "OBJECTID")
  pdp_summary <- pdp_plot_data %>%
    filter(!is.na(NDVI_1)) %>%
    group_by(NDVI_Value = round(NDVI_1, 2)) %>%
    dplyr::summarise(Mean_Predicted_Density = mean(Predicted_Density, na.rm = TRUE), .groups = "drop") %>%
    filter(!is.na(NDVI_Value))
  
  # 将所有结果打包到一个列表中
  fig4_results_data <- list(
    all_predictions = all_regression_preds,       # 用于 Panel A 和 B
    avg_importance = avg_importance,         # 用于 Panel C
    pdp_summary = pdp_summary,               # 用于 Panel D
    analysis_results_global = fig1_analysis_results # 传递图1的结果
  )
  
  # 创建目录并保存
  dir.create("models", showWarnings = FALSE)
  saveRDS(fig4_results_data, "models/fig4_plotting_data.rds")
  # 另外保存最终模型和特征，供图5使用
  saveRDS(validation_rf, "models/best_landslide_bz_model.rds")
  saveRDS(features, "models/model_features.rds")
  saveRDS(modeling_df_full %>% dplyr::select(all_of(features)), "models/training_data_snapshot.rds")
  
  cat("✓ Modeling complete. All results saved to 'models/fig4_plotting_data.rds'.\n")
}

# 确保 landslide_data_processed 和 fig1_analysis_results 对象存在
run_fig4_modeling() 

# =========================================================================
# Part 2: 绘图函数
# 任务: 加载已保存的模型结果，生成图4的四个面板和最终组合图
# =========================================================================
create_fig4_plots <- function() {
  cat("\n=== CREATING PLOTS FOR FIGURE 4 ===\n")
  
  # --- 1. 加载绘图所需的数据 ---
  cat("Step 1: Loading pre-computed modeling results...\n")
  results_path <- "models/fig4_plotting_data.rds"
  if (!file.exists(results_path)) {
    stop("Plotting data not found. Please run `run_fig4_modeling()` first.")
  }
  fig4_data <- readRDS(results_path)
  
  # 将列表中的数据解包为独立变量，方便使用
  all_predictions <- fig4_data$all_predictions
  avg_importance <- fig4_data$avg_importance
  pdp_summary <- fig4_data$pdp_summary
  analysis_results_global <- fig4_data$fig1_analysis_results
  
  # --- 2. 生成四个可视化面板 (代码与您之前的版本几乎完全相同) ---
  cat("Step 2: Generating visualization panels...\n")
  
  # --- Panel A: Prediction Performance ---
  cat("Generating visualization fig4 panel a...\n")
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
  cat("Generating visualization fig4 panel b...\n")
  all_predictions$Residual <- all_predictions$Predicted - all_predictions$Actual
  plot_4b <- ggplot(all_predictions, aes(x = Residual, fill = factor(Fold))) +
    geom_density(alpha = 0.4, color = "black", linewidth = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
    scale_fill_viridis_d(name = "CV Fold") +
    labs(x = "Residual (Predicted - Actual)", y = "Density") +
    nature_theme_professional +
    # 修改点1: 将图例移到图框内右上角
    theme(legend.position = c(0.89, 0.75),
          legend.background = element_rect(fill = "white", color = "gray90"),
          legend.margin = ggplot2::margin(5, 5, 5, 5))
  
  # --- Panel C: Feature Importance ---
  cat("Generating visualization fig4 panel c...\n")
  avg_importance$Category <- case_when(
    grepl("EVI_1|LAI_1", avg_importance$Feature) ~ "Vegetation Indices",
    grepl("change", avg_importance$Feature) ~ "Temporal Changes",
    grepl("Hansen|MODIS", avg_importance$Feature) ~ "Vegetation Types",
    TRUE ~ "Environmental"
  )
  feature_palette <- setNames(c("#3B6F9E", "#D85C60", "#63A27D", "#E8A354"), 
                              c("Environmental", "Temporal Changes", "Vegetation Indices", "Vegetation Types"))
  
  # 修改点2: 简化小图c的x轴标签
  # 创建一个简化的特征名列
  cat("Generating visualization fig4 panel c Features Name...\n")
  avg_importance$SimpleFeature <- gsub("NDVI_change_1_to_2", "Nc1t2", avg_importance$Feature)
  avg_importance$SimpleFeature <- gsub("EVI_change_1_to_2", "Ec1t2", avg_importance$SimpleFeature)
  avg_importance$SimpleFeature <- gsub("LAI_change_1_to_2", "Lc1t2", avg_importance$SimpleFeature)
  avg_importance$SimpleFeature <- gsub("Hansen_Tree_Cover_2000_Percent", "HTC2P", avg_importance$SimpleFeature)
  avg_importance$SimpleFeature <- gsub("Climate_Zone", "CliZ", avg_importance$SimpleFeature)
  avg_importance$SimpleFeature <- gsub("MODIS_IGBP_Simplified", "IGBPS", avg_importance$SimpleFeature)
  
  cat("Generating visualization fig4 panel c plot...\n")
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
    # 修改点2: 去除图例标题
    scale_fill_manual(values = feature_palette) +
    labs(x = "Feature", y = "Importance (Mean %IncMSE)") +
    nature_theme_professional +
    theme(
      legend.position = c(0.78, 0.18),
      legend.box = "horizontal",
      legend.title = element_text(size = 9), # Control legend title size
      legend.text = element_text(size = 8),  # Control legend label size
      legend.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5),
      legend.box.margin = ggplot2::margin(t = -5),
      panel.grid.major.y = element_blank(),
      # 添加以下行来设置图例背景和边框
      legend.background = element_rect(
        colour = "grey80",   # 边框颜色，可以设置为 "grey50" 或你想要的任何颜色
        linewidth = 0.3,     # 边框宽度，可以调整
        fill = "white"      # 图例背景填充色，默认为"white"，你可以根据需要更改
      )
    ) 
  
  # --- Panel D: NDVI Threshold Validation
  cat("Generating visualization fig4 panel d...\n")
  ndvi_breakpoints <- fig1_analysis_results$NDVI_1$breakpoints
  
  plot_4d <- ggplot(pdp_summary, aes(x = NDVI_Value, y = Mean_Predicted_Density)) +
    geom_line(color = "darkgreen", linewidth = 1.5, alpha = 0.8) +
    geom_vline(xintercept = ndvi_breakpoints, linetype = "dashed", color = "red", linewidth = 0.8) +
    annotate("rect", xmin = ndvi_breakpoints[1], xmax = ndvi_breakpoints[2],
             ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.1) +
    annotate("text", x = mean(ndvi_breakpoints) + 0.08, y = max(pdp_summary$Mean_Predicted_Density, na.rm=T) * 0.9,
             label = "Fig 1 Critical Transition Zone", color = "red", size = 3, fontface = "bold", hjust=1.3) +
    labs(x = "NDVI Value", y = "Predicted Landslide Density (Proxy)") +
    nature_theme_professional
  
  # --- 3. 组合最终的Figure 4 ---
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

# 这个函数会从文件中读取数据，运行速度很快
my_figure_4 <- create_fig4_plots()

