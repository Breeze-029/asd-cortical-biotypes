# ==============================================================================
# 图B 高清重出（2026-10-12）：为 Fig4 新排版（A 顶部整行 / B 左下 4×2 / C 右下）
# 用途：
#   1) b_grid_HD.png   2行×4列 脑图网格（无图例），供 Python 端按行优先重排为 4×2
#   2) b_legend_HD.png 同一图带底部色标，供裁出色标条
# 参数与 Fig3_03_b_brainPAD_region_maps.R 完全一致（r 范围 ±sym_limit，
# RdBu 风格 low=#0072B2 mid=white high=#D55E00，grey78 分割线）。
# ==============================================================================
suppressMessages({
  library(ggseg); library(ggplot2); library(dplyr); library(readxl); library(scales)
})

# 2026-10-15：数据源由未校正 PAD 的 ASD_All_ROI_PAD_Correlation_Wscore.xlsx
#   改为 GPR 校正版（tmpRunFig301GPR.m 输出），以与面板标题 "corrected Brain-PAD"
#   及正文数值（左楔前叶 r=−0.422、61 of 68 区域 FDR 显著）一致。
# ---- Repository paths (no edit needed if you keep this folder layout) ----
script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir  <- if (length(script_file)) dirname(normalizePath(sub("^--file=", "", script_file[1]))) else normalizePath(".")
repo_root   <- dirname(script_dir)
data_dir    <- file.path(repo_root, "data")
out_dir     <- file.path(repo_root, "output")
dir.create(out_dir, showWarnings = FALSE)

file_path <- file.path(data_dir, "GPR_ASD_All_ROI_PAD_Correlation_Wscore.xlsx")
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

sym_limit <- max(0.1, max(abs(plot_data$display_val), na.rm = TRUE))

base_theme <- theme_void() +
  theme(
    text = element_text(family = "sans"),
    legend.title = element_text(size = 6, face = "bold"),
    legend.text = element_text(size = 5),
    legend.position = "bottom",
    legend.margin = margin(t = 0.2, unit = "cm"),
    plot.margin = unit(c(0.1, 0.2, 0.2, 0.2), "cm"),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  )

make_plot <- function(show_legend) {
  ggplot(plot_data) +
    geom_brain(
      atlas = dk(),
      position = position_brain(hemi ~ view),
      aes(fill = display_val),
      colour = "grey78", size = 0.1
    ) +
    scale_fill_gradient2(
      low = "#0072B2", mid = "white", high = "#D55E00",
      midpoint = 0, limits = c(-sym_limit, sym_limit),
      oob = squish, na.value = "grey95", name = "r value"
    ) +
    # 色标长度统一为脑图网格宽的 ~40%（与 Fig2 Subtype 脑图比例一致），厚 2.2 mm
    guides(fill = guide_colourbar(barwidth = unit(32, "mm"),
                                  barheight = unit(2.2, "mm"))) +
    base_theme +
    (if (show_legend) NULL else theme(legend.position = "none"))
}

# 高分辨率输出（9.5cm × 7.5cm，900 dpi ≈ 3366 × 2657）
ggsave(file.path(out_dir, "b_grid_HD.png"),   make_plot(FALSE), width = 9.5, height = 7.5,
       units = "cm", dpi = 900, bg = "white")
ggsave(file.path(out_dir, "b_legend_HD.png"), make_plot(TRUE),  width = 9.5, height = 7.5,
       units = "cm", dpi = 900, bg = "white")
cat("sym_limit =", sym_limit, "\n")
cat("saved", file.path(out_dir, "b_grid_HD.png"), ",", file.path(out_dir, "b_legend_HD.png"), "\n")
