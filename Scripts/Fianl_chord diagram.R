rm(list = ls())
setwd("E:\\ymkim\\experiment file\\Gut-Brain axis\\종근당 바이오\\3차\\JONG_3RD_SAMPLE\\MZMINE\\Known data(논문용)\\Correlation\\Kynurenic acid with others")

library(circlize)
library(reshape2)
library(tidyr)
library(dplyr)

Sys.setlocale(category = "LC_ALL", locale = "Greek")

# 데이터 불러오기
mat <- as.data.frame(readxl::read_excel("[Result] Total_metabolite_correlation.xlsx", sheet = "Chord diagram_kynurenic acid_cS"))

nm <- unique(c(unlist(mat["M1.TYPE1"]), unlist(mat["M2.TYPE2"])))
group <- structure(gsub("_\\d*", "", nm), names = nm)

grid.col <- group
grid.col <- replace(grid.col, grid.col == "Cecum", "#A0FFA0")
grid.col <- replace(grid.col, grid.col == "Hp", "skyblue")
grid.col <- replace(grid.col, grid.col == "Serum", "#F9CCF0")
# sector 이름을 metabolite 이름으로 사용
re_nm <- unique(c(as.character(mat$Metabolite1), as.character(mat$Metabolite2)))
names(grid.col) <- re_nm

re_mat <- mat[, c(3,7,9)]
col_fun <- colorRamp2(c(-1, 0, 1), c("#55a0fb", "white", "#FF8080"), transparency = 0.3)
re_mat$color <- col_fun(re_mat$coefficient)

re_group <- structure(gsub("_\\d*", "", nm), names = re_nm)
group_order <- c("Serum", "Hp", "Cecum")
re_group_02 <- re_group[order(match(re_group, group_order))]

grid.col_02 <- re_group_02
grid.col_02 <- replace(grid.col_02, grid.col_02 == "Cecum", "#008000")
grid.col_02 <- replace(grid.col_02, grid.col_02 == "Serum", "#9F044D")
grid.col_02 <- replace(grid.col_02, grid.col_02 == "Hp", "#0351AF")
#grid.col_02 <- replace(grid.col_02, grid.col_02 == "Feces", "#0000AF")

for.highlight <- as.data.frame(grid.col_02)
for.highlight$highlight <- rownames(for.highlight)

# metabolite 빈도수 계산
metabolite_all <- c(as.character(mat$Metabolite1), as.character(mat$Metabolite2))


# PDF 저장
pdf(file = "02_hard_version_CS.pdf", width = 15, height = 15)
circos.clear()

circos.par(canvas.ylim = c(-2.7, 2.7), 
           track.margin = c(0.01, 0.01),
           track.height = 0.05,
           cell.padding = c(1, 1, mm_h(10), 1))

# 섹터 이름을 기반으로 투명한 색상 벡터 생성
grid.col.transparent <- setNames(rep("transparent", length(re_nm)), re_nm)
# chord diagram 그리기 (내부 섹터 및 라벨 표시)
chordDiagram(
  re_mat,
  col = re_mat[, 4],
  group = re_group_02,
  annotationTrack = "grid",
  annotationTrackHeight = 0.05,
  grid.col = grid.col.transparent ,
  preAllocateTracks = list(track.height = mm_h(4), track.margin = c(mm_h(0.5), 0)),
  small.gap = 1.5,
  big.gap = 7
)

# histogram 트랙 추가 (track.index = 2, y축 범위 0~3)
circos.trackPlotRegion(
  track.index = 2,
  ylim = c(0,3),
  col = grid.col[metabolite]
)

# sector label track (track.index = 1)
circos.track(track.index = 1, panel.fun = function(x, y) {
  # 원래 섹터 이름에서 접두어 제거 (정규표현식 사용)
  label <- sub("^(Hp_|Cecum_|Serum_)", "", CELL_META$sector.index)
  circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 1.5, 
              label, 
              facing = "clockwise", 
              niceFacing = FALSE,
              adj = c(0, 0.5), 
              cex = 0.8)
}, bg.border = NA)

# highlight 예시 (필요시)
highlight.sector(rownames(for.highlight[for.highlight$grid.col == "#008000",]),
                 track.index = 1,
                 col = "#008000",
                 text = "Cecum", cex = 0.9, text.col = "white", 
                 facing = "bending.inside", niceFacing = TRUE)

highlight.sector(rownames(for.highlight[for.highlight$grid.col == "#9F044D",]),
                 track.index = 1,
                 col = "#9F044D",
                 text = "Serum", cex = 0.9, text.col = "white", 
                 facing = "bending.inside", niceFacing = TRUE)

highlight.sector(rownames(for.highlight[for.highlight$grid.col == "#0351AF",]),
                 track.index = 1,
                 col = "#0351AF",
                 text = "Hp", cex = 0.9, text.col = "white", 
                 facing = "bending.inside", niceFacing = TRUE)



circos.clear()

dev.off()

