library(readr)
library(ggplot2)
library(dplyr)
library(rtracklayer)
library(reshape2)
library(VennDiagram)
library(RColorBrewer)
library(igraph)


data_dir <- file.path("data", "subgenome")
out_dir <- file.path("analysis", "subgenome")

HCoV_chrom_size <- 27317
HCoV_sample <- c("WT", "SL2")
HCoV_nanopore_bw_WT <- import(file.path(data_dir, "HCoV", "bw", "Nanopore.WT.Forward.bw"), "bw")
HCoV_nanopore_bw_SL2 <- import(file.path(data_dir, "HCoV", "bw", "Nanopore.SL2.Forward.bw"), "bw")

################################################
## Figure S8C
position_df <- data.frame(Position=1:HCoV_chrom_size)
signal_WT <- data.frame(Position=HCoV_nanopore_bw_WT@ranges@start, Signal=HCoV_nanopore_bw_WT$score)
signal_WT <- left_join(position_df, signal_WT)
signal_WT$Sample <- "WT"
signal_SL2 <- data.frame(Position=HCoV_nanopore_bw_SL2@ranges@start, Signal=HCoV_nanopore_bw_SL2$score)
signal_SL2 <- left_join(position_df, signal_SL2)
signal_SL2$Sample <- "SL2"
merge_signal_df <- rbind(signal_WT, signal_SL2)
merge_signal_df$Signal[is.na(merge_signal_df$Signal)] <- 0
p <- ggplot(merge_signal_df, aes(x=Position, y=log10(Signal+1), color=Sample)) +
  geom_line(size=0.1) +
  theme_bw() +
  scale_y_continuous(breaks = log10(c(0+1, 1+1, 10+1, 100+1, 1000+1, 10000+1)), labels = c("0", "1E0", "1E1", "1E2", "1E3", "1E4")) +
  labs(y="Norm. signal (HCoV)") +
  theme(text = element_text(family="ArialMT", color = "black", size = 6),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.background = element_blank(),
        legend.title = element_text(family="ArialMT", color = "black", size = 6),
        legend.key.size = unit(4, "mm"),
        panel.grid = element_blank())
ggsave(file.path(oION, "Nanopore.singal.profile.pdf"), p, width = 7, height = 3.5, units = "cm")
