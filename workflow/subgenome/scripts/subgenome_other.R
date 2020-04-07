library(readr)
library(ggplot2)
library(dplyr)


data_dir <- file.path("data", "subgenome")
out_dir <- file.path("analysis", "subgenome")
if (!dir.exists(out_dir)) dir.create(out_dir)

polyA_tail_rep1 <- read_delim("data/polyA_tail/Nanopore_rep1.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
polyA_tail_rep2 <- read_delim("data/polyA_tail/Nanopore_rep2.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

nanopore_bed_rep1 <- read_delim(file.path(data_dir, "nanopore", "data", "WIV04.rep1.bed"), 
                                "\t", escape_double = FALSE, col_names = FALSE,
                                col_types = cols(X10 = col_integer(), 
                                                 X11 = col_character(), X12 = col_character(), 
                                                 X7 = col_integer(), X8 = col_integer(), 
                                                 X9 = col_character()), trim_ws = TRUE)
nanopore_bed_rep2 <- read_delim(file.path(data_dir, "nanopore", "data", "WIV04.rep2.bed"), 
                                "\t", escape_double = FALSE, col_names = FALSE,
                                col_types = cols(X10 = col_integer(), 
                                                 X11 = col_character(), X12 = col_character(), 
                                                 X7 = col_integer(), X8 = col_integer(), 
                                                 X9 = col_character()), trim_ws = TRUE)

mean(nanopore_bed_rep1$X4 %in% polyA_tail_rep1$ReadID)
mean(nanopore_bed_rep2$X4 %in% polyA_tail_rep2$ReadID)

polyA_tail_rep1 <- polyA_tail_rep1[polyA_tail_rep1$ReadID %in% nanopore_bed_rep1$X4,]
polyA_tail_rep2 <- polyA_tail_rep2[polyA_tail_rep2$ReadID %in% nanopore_bed_rep2$X4,]
polyA_tail <- rbind(polyA_tail_rep1, polyA_tail_rep2)
p <- ggplot(polyA_tail, aes(x=Length)) +
  geom_histogram(fill="black", bins=100) +
  theme_bw() +
  labs(x="polyA-tail length", y="#Reads") +
  lims(x=c(0, 200)) +
  theme(text = element_text(family="ArialMT", size=6),
        title = element_text(family="ArialMT", size=6),
        axis.text = element_text(color = "black"),
        panel.grid = element_blank()
  )
ggsave(file.path(out_dir, "Nanopore", "polyA.pdf"), p, width = 12, height = 6, units = "cm")

Pair_info_all <- read_delim("analysis/subgenome/seq_pair/Pair.info.all.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
plot_pair <- function(df){
  DL_AL <- df[df$PairType=="DL-AL",]
  DR_AR <- df[df$PairType=="DR-AR",]
  p <- ggplot() +
    geom_curve(DR_AR, mapping = aes(x=Start, y=0, xend=End, yend=0, alpha=log10(ReadNum)), curvature=0.5, color="red", size=0.1) +
    geom_curve(DL_AL, mapping = aes(x=Start, y=0, xend=End, yend=0, alpha=log10(ReadNum)), curvature=-0.5, color="blue", size=0.1) +
    theme_bw() +
    labs(x="Position") +
    theme(text = element_text(family="ArialMT", size=6),
          title = element_text(family="ArialMT", size=6),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),
          axis.text = element_text(color = "black"),
          legend.key.size = unit(4, "mm"),
          # legend.title = element_blank(),
          legend.background = element_blank(),
          legend.text = element_text(family="ArialMT", size=6),
          panel.grid = element_blank()
    )
  return(p)
}

p <- plot_pair(Pair_info_all)
ggsave(file.path(out_dir, "seq_pair", "JS.gap_type.curve.all.pdf"), p, width = 12, height = 6, units = "cm")

consistent_df <- Pair_info_all[Pair_info_all$IsConsistent,]
p <- plot_pair(consistent_df)
ggsave(file.path(out_dir, "seq_pair", "JS.gap_type.curve.consistent.pdf"), p, width = 12, height = 6, units = "cm")

