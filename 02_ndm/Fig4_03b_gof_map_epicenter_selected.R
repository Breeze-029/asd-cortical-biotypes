# ============================================================================
# GOF map — 震源区高亮版（GOF map epicenter selected）
#   仅给 NDM 识别出的震源 ROI 上色，其余脑区为白色
# 输入：AllROI_FullGOF_Volume.xlsx（列：ROI, SUB1, SUB2；行 = ROI 1–68，奇左偶右）
# 输出：Brain_Maps_GOF/GOF_map_epicenter_selected_SUB1.pdf(.png)、_SUB2.pdf(.png)
# 规范：与 GOF_map.R 完全一致（89.5 × 70.7 mm、Arial、RdBu 蓝-白-红、矢量 PDF + 600 dpi）
# 修复：原实现按“区域名”匹配导致左右半球同时高亮；
#       现按 (ROI 编号 → 半球) 精确匹配：ROI 奇 = 左半球，偶 = 右半球
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

# NDM 识别出的震源区（ROI 编号）：按正文规则取 GOF 前 15%（按排名，round(0.15*n) = 10 个区域）
n_top <- round(0.15 * nrow(data))
epi_list <- lapply(c("SUB1", "SUB2"), function(s) {
  v <- data[[s]]
  data$ROI[order(-v)][seq_len(n_top)]
})
names(epi_list) <- c("SUB1", "SUB2")
message("TOP15% n = ", n_top)
message("SUB1 GOF 前15% 区域(ROI): ", paste(sort(epi_list$SUB1), collapse = ", "))
message("SUB2 GOF 前15% 区域(ROI): ", paste(sort(epi_list$SUB2), collapse = ", "))
message("SUB1 GOF 前15% 区域(ROI): ", paste(epi_list$SUB1, collapse = ", "))
message("SUB2 GOF 前15% 区域(ROI): ", paste(epi_list$SUB2, collapse = ", "))

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
mm_w <- 8.95   # cm（89.5 mm）
mm_h <- 7.07   # cm（70.7 mm）

# 与 GOF_map.R 使用同一色阶，保证两套图可直接比较
all_vals  <- c(data$SUB1, data$SUB2)
limit_val <- max(abs(range(all_vals, na.rm = TRUE))) * 1.02

# 关键修复：按 ROI 编号（含半球信息）精确匹配，只高亮对应半球
prepare_brain_data <- function(gof_col, epi_rois) {
  data.frame(ROI = 1:68, gof = gof_col) %>%
    mutate(hemi       = if_else(ROI %% 2 == 1, "left", "right"),
           region_idx = ceiling(ROI / 2),
           region     = dk_region_names[region_idx],
           is_epi     = ROI %in% epi_rois,
           gof_display = ifelse(is_epi, gof, NA_real_)) %>%
    select(region, hemi, gof, gof_display, is_epi)
}

for (sub in c("SUB1", "SUB2")) {
  message(sprintf("正在绘制震源高亮 GOF 图: %s", sub))
  plot_df <- prepare_brain_data(data[[sub]], epi_list[[sub]])

  # 打印被高亮的半球-区域及 GOF 值，便于核对（左右半球分开列出）
  sel_tbl <- plot_df %>% filter(is_epi) %>% select(region, hemi, gof)
  message(sprintf("[%s] 高亮 %d 个（半球, 区域）组合：", sub, nrow(sel_tbl)))
  print(as.data.frame(sel_tbl))

  p <- ggplot(plot_df) +
    geom_brain(atlas = dk(), position = position_brain(hemi ~ view),
               aes(fill = gof_display),
               colour = "grey78", size = 0.1) +
    scale_fill_gradient2(
      # RdBu 蓝-白-红，与图1（Subtype mean cortical volume W-score 图）色标一致
      low = "#2166ac", mid = "white", high = "#b2182b",
      midpoint = 0,
      limits = c(-limit_val, limit_val),
      oob = squish,
      na.value = "grey95"     # 非震源区：浅灰底（与图1 一致），避免纯白背景下孤立色块观感偏浅
    ) +
    # 色标长度统一为脑图网格宽的 ~40%（与 Fig2 Subtype 脑图比例一致），厚 2.2 mm
    guides(fill = guide_colourbar(barwidth = unit(31, "mm"),
                                  barheight = unit(2.33, "mm"))) +
    theme_void() +
    theme(
      text = element_text(family = "sans"),
      legend.position = "bottom",
      legend.title = element_blank(),
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

  pdf_file <- file.path(output_dir, paste0("GOF_map_epicenter_selected_", sub, ".pdf"))
  png_file <- file.path(output_dir, paste0("GOF_map_epicenter_selected_", sub, ".png"))
  ggsave(filename = pdf_file, plot = p, width = mm_w, height = mm_h, units = "cm", device = pdf)
  ggsave(filename = png_file, plot = p, width = mm_w, height = mm_h, units = "cm", dpi = 600)
  message(sprintf("已输出: %s", pdf_file))
}

message("GOF map (epicenter selected) 完成：2 张，各 89.5 × 70.7 mm")
