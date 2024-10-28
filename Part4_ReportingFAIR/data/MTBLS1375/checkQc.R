library(tidyverse)
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('PCAtools')
# https://bioconductor.org/packages/devel/bioc/vignettes/PCAtools/inst/doc/PCAtools.html
library('PCAtools')
skyl_col = cols(
  OriginalTransitionName = col_character(),
  AnalyticalID = col_character(),
  PrecursorMz = col_double(),
  ProductMz = col_double(),
  `Fragment Ion` = col_character(),
  RT = col_double(),
  Area = col_double(),
  `Product Ion Formula` = col_character(),
  `Precursor Note` = col_character(),
  `Peptide Note` = col_character(),
  `Transition Replicate Note` = col_character()
)
sampleInfo <- read_tsv("s_MTBLS1375.txt")
skyl <-
  read_csv("human_plasma_export.csv",
           na = c("NA", "#N/A", ""),
           col_types = skyl_col) %>%
  left_join(sampleInfo, by=c("AnalyticalID"="Source Name")) %>%
  
  
  mutate("OriginalTransitionName" = str_replace(.$"OriginalTransitionName", " ([abc])$", "\\1")) %>%
  mutate(
    "OriginalTransitionName" = str_replace(
      .$"OriginalTransitionName",
      "^([A-Zderx0-9 :;\\-_\\/\\(\\)\\[\\]]+)([abc]?)$",
      "\\1"
    ),
    "isoid" = str_replace(
      .$"OriginalTransitionName",
      "^([A-Zderx0-9 :;\\-_\\/\\(\\)\\[\\]]+)([abc]?)$",
      "\\2"
    )
  ) %>%
  mutate("OriginalTransitionName" = str_replace(.$"OriginalTransitionName", "^(.*)[ ]?$", "\\1")) %>%
  select("OriginalTransitionName", "isoid", everything()) %>% group_by(OriginalTransitionName, AnalyticalID, PrecursorMz, ProductMz, `Fragment Ion`, RT) %>%
  summarise(Area = sum(Area))

ggplot(skyl %>% filter(OriginalTransitionName=="CE 14:0"), aes(x=fct_reorder(AnalyticalID, Area), y=Area, fill=`Characteristics[Sample type]`)) + geom_boxplot() + facet_grid(.~`Characteristics[Sample type]`, drop = FALSE, scales = "free") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

df <- skyl %>% pivot_wider(id_cols="AnalyticalID", names_from="OriginalTransitionName", values_from = "Area", values_fn = sum) %>% ungroup()
mim <- as.matrix(df %>% select(-AnalyticalID)) %>% replace_na(0)
rownames(mim) <- df$AnalyticalID
library(pheatmap)
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
mim_nm <- t(apply(t(mim), 1, cal_z_score))
mim_nm <- na.omit(mim_nm)

pdf("plasma-heatmap-raw-areas.pdf", width = 20, height = 44)
#p <- pheatmap(mim_nm, cutree_cols = 8, cutree_rows = 23, fontsize_row = 8, fontsize_col = 10, width = 8, height = 12)
p <- pheatmap(mim_nm, cutree_cols = 27, cutree_rows = 30, fontsize_row = 8, fontsize_col = 10, width = 12, height = 20)
grid::grid.draw(p)
dev.off()

pca <- princomp(min_nm)