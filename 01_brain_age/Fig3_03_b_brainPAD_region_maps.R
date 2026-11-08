# ==============================================================================
# 图B 全脑Brain-PAD与形态偏差关联图（去总标题精简版 | 统一AI加字）
# ==============================================================================
library(ggseg)
library(ggplot2)
library(dplyr)
library(readxl)
library(scales)

# ---- Repository paths (no edit needed if you keep this folder layout) ----
script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir  <- if (length(script_file)) dirname(normalizePath(sub("^--file=", "", script_file[1]))) else normalizePath(".")
repo_root   <- dirname(script_dir)
data_dir    <- file.path(repo_root, "data")
out_dir     <- file.path(repo_root, "output")
dir.create(out_dir, showWarnings = FALSE)

file_path <- file.path(data_dir, "ASD_All_ROI_PAD_Correlation_Wscore.xlsx")

if (!file.exists(file_path)) {
  stop("【错误】找不到相关性结果文件，请检查路径！")
}

my_data <- read_xlsx(file_path)

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

plot_data <- my_data %>%
  mutate(
    hemi = if_else(row_number() %% 2 == 1, "left", "right"),
    region_idx = ceiling(row_number() / 2),
    region = dk_region_names[region_idx],
    display_val = Correlation_R
  )

max_abs <- max(abs(plot_data$display_val), na.rm = TRUE)
sym_limit <- max(0.1, max_abs)

brain_plot <- ggplot(plot_data) +
  geom_brain(
    atlas = dk(), 
    position = position_brain(hemi ~ view), 
    aes(fill = display_val),
    colour = "grey78", size = 0.1
  ) + 
  scale_fill_gradient2(
    low = "#0072B2",
    mid = "white",
    high = "#D55E00",
    midpoint = 0,
    limits = c(-sym_limit, sym_limit),
    oob = squish,
    na.value = "grey95",    
    name = "r value"
  ) +
  theme_void() + 
  theme(
    # 兼容Mac系统字体，避免报错
    text = element_text(family = "sans"),
    # ==========【修改1：移除ggtitle总标题，统一到AI中添加】==========
    # 图例保留，色标是图形核心元素，AI无需重做
    legend.title = element_text(size = 6, face = "bold"),
    legend.text = element_text(size = 5),
    legend.key.width = unit(1.2, "cm"),
    legend.key.height = unit(2.2, "mm"),
    legend.position = "bottom",
    legend.margin = margin(t = 0.2, unit = "cm"),
    # ==========【修改2：大幅收紧上下边距，消除多余空白】==========
    # 上、右、下、左；顶部几乎不留空，底部仅留图例空间
    plot.margin = unit(c(0.1, 0.2, 0.2, 0.2), "cm"),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  )

# 导出设置：宽9.5cm，高7.5cm，刚好容纳脑图+底部图例，无多余空白
output_dir <- out_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

pdf_file <- file.path(output_dir, "Figure_B_BrainPAD_Regional_Correlations_noTitle.pdf")
png_file <- file.path(output_dir, "Figure_B_BrainPAD_Regional_Correlations_noTitle.png")

ggsave(
  filename = pdf_file,
  plot = brain_plot, 
  width = 9.5, height = 7.5, units = "cm", 
  dpi = 300, device = "pdf"
)
ggsave(
  filename = png_file,
  plot = brain_plot, 
  width = 9.5, height = 7.5, units = "cm", 
  dpi = 300
)

print(brain_plot)
cat("✅ 图B导出完成，尺寸：9.5cm × 7.5cm（无总标题，仅脑图+色标图例）\n")
cat("颜色范围：", -sym_limit, " 至 ", sym_limit, "\n")
