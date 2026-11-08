# ============================================================================
# 全脑 GOF（NDM 拟合优度）皮层分布图 — Desikan–Killiany 68 区
# 输入：AllROI_FullGOF_Volume.xlsx（列：ROI, SUB1, SUB2；行 = ROI 1–68，奇左偶右）
# 输出：Brain_Maps_GOF/GOF_map_SUB1.pdf(.png)、GOF_map_SUB2.pdf(.png)
# 规范：Arial(sans)；单张 89.5 × 70.7 mm（183 mm 宽、行内两亚型并排、共三行布局）；
#       矢量 PDF + 600 dpi PNG；色标与图1（Subtype mean cortical volume W-score 图）
#       一致（RdBu 蓝-白-红）；无色条标题、无标题/说明文字（留给 AI 与图注）
# ============================================================================
suppressMessages({
  library(ggseg); library(ggplot2); library(dplyr); library(readxl); library(scales)
})

# ---- Repository paths (no edit needed if you keep this folder layout) ----
script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir  <- if (length(script_file)) dirname(normalizePath(sub("^--file=", "", script_file[1]))) else normalizePath(".")
repo_root   <- dirname(script_dir)
data_dir    <- file.path(repo_root, "data")
out_dir     <- file.path(repo_root, "output")
dir.create(out_dir, showWarnings = FALSE)

file_path  <- file.path(data_dir, "Fig4_allROI_GOF_volume.xlsx")
output_dir <- file.path(repo_root, "output", "Brain_Maps_GOF")
if (!file.exists(file_path)) stop("找不到 Excel 文件: ", file_path)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

data <- read_xlsx(file_path, sheet = 1)
names(data)[1:3] <- c("ROI", "SUB1", "SUB2")

# DK 34 区标准名称（ROI 1–68，奇 = 左半球，偶 = 右半球）
dk_region_names <- c(
  "precentral", "postcentral", "paracentral", "superior frontal",
  "rostral middle frontal", "caudal middle frontal", "frontal pole",
  "lateral orbitofrontal", "medial orbitofrontal", "pars opercularis",
  "pars orbitalis", "pars triangularis", "insula",
  "rostral anterior cingulate", "caudal anterior cingulate", "posterior cingulate",
  "isthmus cingulate", "parahippocampal", "superior parietal", "inferior parietal",
  "cuneus", "entorhinal", "fusiform", "lingual",
  "pericalcarine", "precuneus", "supramarginal", "superior temporal",
  "middle temporal", "inferior temporal", "temporal pole",
  "banks of superior temporal sulcus", "transverse temporal", "lateral occipital"
)

# 拼版尺寸：183 mm 宽、共三行（每组一行，行内两亚型并排，中缝 4 mm）
mm_w <- 8.95   # cm（89.5 mm = (183 − 4)/2）
mm_h <- 7.07   # cm（70.7 mm，保持 DK 脑图 1.266:1 比例）

# 两个亚型共用同一对称色阶（以 0 为中心），保证两张图可直接相互比较
all_vals  <- c(data$SUB1, data$SUB2)
limit_val <- max(abs(range(all_vals, na.rm = TRUE))) * 1.02

prepare_brain_data <- function(gof_col) {
  data.frame(ROI = 1:68, gof = gof_col) %>%
    mutate(hemi = if_else(ROI %% 2 == 1, "left", "right"),
           region_idx = ceiling(ROI / 2),
           region = dk_region_names[region_idx]) %>%
    select(region, hemi, gof)
}

for (sub in c("SUB1", "SUB2")) {
  message(sprintf("正在绘制全脑 GOF 图: %s", sub))
  plot_df <- prepare_brain_data(data[[sub]])

  p <- ggplot(plot_df) +
    geom_brain(atlas = dk(), position = position_brain(hemi ~ view), aes(fill = gof),
               colour = "grey78", size = 0.1) +
    scale_fill_gradient2(
      # RdBu 蓝-白-红，与图1（Subtype mean cortical volume W-score 图）色标一致
      low = "#2166ac",    # 低/负：蓝色
      mid = "white",      # 0：白色
      high = "#b2182b",   # 高/正：红色
      midpoint = 0,
      limits = c(-limit_val, limit_val),
      oob = squish,
      na.value = "grey95"      # 无数据区：与图1（Subtype W-score 图）一致
    ) +
    # 色标长度统一为脑图网格宽的 ~40%（与 Fig2 Subtype 脑图比例一致），厚 2.2 mm
    guides(fill = guide_colourbar(barwidth = unit(31, "mm"),
                                  barheight = unit(2.33, "mm"))) +
    theme_void() +
    theme(
      text = element_text(family = "sans"),
      legend.position = "bottom",
      legend.title = element_blank(),          # 色条标题留给 AI/图注统一添加
      legend.text = element_text(size = 5),
      legend.margin = margin(t = 0.15, unit = "cm"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      plot.margin = unit(c(0.1, 0.15, 0.15, 0.15), "cm"),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.background = element_rect(fill = "white", color = NA)
    )

  pdf_file <- file.path(output_dir, paste0("GOF_map_", sub, ".pdf"))
  png_file <- file.path(output_dir, paste0("GOF_map_", sub, ".png"))
  # 使用基础 pdf() 设备（本机 arm64 R 无 X11/cairo，cairo_pdf 不可用）
  ggsave(filename = pdf_file, plot = p, width = mm_w, height = mm_h, units = "cm", device = pdf)
  ggsave(filename = png_file, plot = p, width = mm_w, height = mm_h, units = "cm", dpi = 600)
  message(sprintf("已输出: %s", pdf_file))
}

message("GOF map 完成：2 张，各 89.5 × 70.7 mm")
