# 1. Identify significant junctions

library(readr)
library(ggplot2)
library(dplyr)
library(rtracklayer)
library(reshape2)
library(VennDiagram)
library(RColorBrewer)
library(igraph)

sample_df <- data.frame(project=c("NarryKim"), sample=c("VeroInf24h"))
data_dir <- file.path("data", "subgenome", "other_data")
out_dir <- file.path("analysis", "subgenome", "other_data")
if (!dir.exists(out_dir)) dir.create(out_dir)
all_consistent_gap <- read_delim("analysis/subgenome/NGS_Nanopore/all_consistent_gap.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)


load_Nanopore_JS <- function(f_bed, f_gap){
  bed <- read_delim(f_bed, "\t", escape_double = FALSE, col_names = FALSE, col_types = cols(X11 = col_character()), trim_ws = TRUE)
  bed <- bed[,c(2,3,4)]
  names(bed) <- c("Start", "End", "ReadName")
  multi_map_name <- unique(bed$ReadName[duplicated(bed$ReadName)])
  bed <- bed[! bed$ReadName %in% multi_map_name,]
  fl_bed <- bed[(bed$Start < fl_cutoff) & (bed$End > (chrom_size-fl_cutoff)),]
  
  gap <- read_delim(f_gap, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  gap <- gap[,c(2,3,4)]
  names(gap) <- c("Start", "End", "ReadName")
  gap <- gap[gap$ReadName %in% fl_bed$ReadName,]
  gap_info <- gap %>% group_by(Start, End) %>% summarise(ReadNum=n())
  gap_info$ReadRatio <- gap_info$ReadNum / nrow(fl_bed)
  
  res_li <- list(
    total_reads=nrow(bed),
    fl_reads=nrow(fl_bed),
    gap=gap_info
  )
  return(res_li)
}


for (sample_indx in 1:nrow(sample_df)) {
  proj <- sample_df$project[sample_indx]
  samp <- sample_df$sample[sample_indx]
  
  ###### 2 Nanopore data
  fl_cutoff <- 45 # 65 - 20, criteria for being full-(sub)genome at both ends
  chrom_size <- 29891
  min_len_cutoff <- 100
  nanopore_bw <- import(file.path(data_dir, "bw", proj, sprintf("Nanopore.%s.Forward.bw", samp)), "bw")
  nanopore_bed <- read_delim(file.path(data_dir, "nanopore", "data", proj, sprintf("WIV04.%s.bed", samp)), 
                                  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  nanopore_fl <- nanopore_bed$X4[(nanopore_bed$X2 < fl_cutoff) & (nanopore_bed$X3 > (chrom_size-fl_cutoff))]
  length(nanopore_fl)
  #' Read Nanopore junctions and normalize the signals
  compute_size_factor <- function(x, bw, window_size = 100){
    return(sum(bw$score[abs(bw@ranges@start-x) < window_size], na.rm = T))
  }
  
  load_nanopore_info <- function(f, fl_li, label, window_size=100){
    gap_df <- read_delim(f, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
    names(gap_df) <- c("Chrom", "Start", "End", "Name", "Score", "Strand")
    gap_df$FL <- gap_df$Name %in% fl_li
    gap_df <- gap_df[gap_df$Strand == "+", ]
    Nanopore_JS <- gap_df %>% group_by(FL, Start, End) %>% summarise(ReadNum=n())
    
    Nanopore_JS$StartFactor <- sapply(Nanopore_JS$Start, compute_size_factor, bw=nanopore_bw, window_size=window_size) + 0.01
    Nanopore_JS$EndFactor <- sapply(Nanopore_JS$End, compute_size_factor, bw=nanopore_bw, window_size=window_size) + 0.01
    Nanopore_JS$Factor <- sqrt(Nanopore_JS$StartFactor * Nanopore_JS$EndFactor)
    Nanopore_JS$ReadNormNum <- Nanopore_JS$ReadNum / Nanopore_JS$Factor
    Nanopore_JS$ReadNormNum <- Nanopore_JS$ReadNormNum / median(Nanopore_JS$ReadNormNum)
    Nanopore_JS$Len <- Nanopore_JS$End - Nanopore_JS$Start
    Nanopore_JS$ReadRatio <- Nanopore_JS$ReadNum / sum(Nanopore_JS$ReadNum)
    Nanopore_JS <- Nanopore_JS %>% group_by(FL) %>% mutate(ReadRatioGroup=ReadNum/sum(ReadNum))
    Nanopore_JS$FL <- factor(Nanopore_JS$FL, levels = c(TRUE, FALSE), labels = c("FL", "non-FL"))
    Nanopore_JS$Source <- label
    return(Nanopore_JS)
  }
  
  nanopore_JS <- load_nanopore_info(file.path(data_dir, "nanopore", "data", proj, sprintf("WIV04.%s.gap.bed", samp)), nanopore_fl, samp)
  js_info <- load_Nanopore_JS(file.path(data_dir, "nanopore", "data", proj, sprintf("WIV04.%s.bed", samp)),file.path(data_dir, "nanopore", "data", proj, sprintf("WIV04.%s.gap.bed", samp)))
  tmp_js_df <- js_info$gap
  tmp_js_df$Sample <- samp
  tmp_js_df <- tmp_js_df[,c("Start", "End", "ReadNum")]
  names(tmp_js_df)[3] <- as.character(samp)
  all_consistent_gap <- left_join(all_consistent_gap, tmp_js_df)
  
  oION <- file.path(out_dir, "Nanopore", proj)
  if (!dir.exists(oION)) dir.create(oION)
  write_tsv(nanopore_JS, file.path(oION,  sprintf("cutoff.%s.sourceData.tsv", samp)))
  
  # Require more than 2 Nanopore reads for valid junctions
  nanopore_read_cutoff <- 2
  nanopore_JS$ReadPass <- nanopore_JS$ReadNum > nanopore_read_cutoff
  
  p <- ggplot(nanopore_JS, aes(x=ReadNum, fill=ReadPass)) +
    geom_bar() +
    facet_grid(vars(FL), vars(Source), scale="free_y") +
    scale_x_continuous(limits = c(-0.1, 20)) +
    scale_fill_manual(values = c("TRUE"="red", "FALSE"="grey70")) +
    labs(x="#Reads", y="#Junctions") +
    theme_bw() +
    theme(text = element_text(family="ArialMT", color = "black", size = 6),
          axis.line = element_line(color = "black"),
          axis.text = element_text(color = "black"),
          axis.ticks = element_line(color = "black"),
          legend.position = "none",
          panel.grid = element_blank())
  ggsave(file.path(oION, sprintf("ReadNum.cutoff.%s.pdf", samp)), p, width = 4, height = 6, units = "cm")
  
  # Compute Nanopore local background cutoff
  nanopore_norm_cutoff <- quantile(nanopore_JS$ReadNormNum[!nanopore_JS$ReadPass], 0.99)  # Top 1%
  nanopore_JS$ReadPass <- factor(nanopore_JS$ReadPass, levels = c(TRUE, FALSE), labels = c("High expr.", "Low expr."))
  p <- ggplot(nanopore_JS, aes(x=ReadNormNum, fill=Source)) +
    geom_histogram(bins=50, alpha=0.3, position="identity") +
    geom_vline(xintercept = nanopore_norm_cutoff, size=0.3, lty="dotted") +
    facet_grid(vars(ReadPass), vars(FL), scales="free_y") +
    scale_x_log10() +
    labs(x="Normalized number", y="#Junctions") +
    theme_bw() +
    theme(text = element_text(family="ArialMT", color = "black", size = 6),
          axis.line = element_line(color = "black"),
          axis.text = element_text(color = "black"),
          axis.ticks = element_line(color = "black"),
          legend.position = "none",
          panel.grid = element_blank())
  ggsave(file.path(oION, sprintf("ReadNormNum.cutoff.%s.pdf", samp)), p, width = 8, height = 6, units = "cm")
  
  # Identify Nanopore junctions passing local background and read# cutoff 
  nanopore_JS$ReadPass <- NULL
  enriched_nanopore_JS <- nanopore_JS[
    ((nanopore_JS$ReadNum>nanopore_read_cutoff) & 
       (nanopore_JS$ReadNormNum>nanopore_norm_cutoff) & 
       (nanopore_JS$Len>min_len_cutoff)), ]
  write_tsv(enriched_nanopore_JS, file.path(oION, sprintf("Nanopore.enriched.%s.sourceData.tsv", samp)))
  
  # Show reproduicity between replicates
  p <- ggplot() +
    geom_point(data=nanopore_JS, mapping=aes(x=Start, y=End), color="grey80", size=0.1, alpha=0.3) +
    geom_point(data=enriched_nanopore_JS, mapping=aes(x=Start, y=End, alpha=log10(ReadNormNum)), size=1, color="red") +
    lims(x=c(0, chrom_size), y=c(0, chrom_size)) +
    facet_wrap(~FL) +
    labs(x="Start position", y="End position", alpha="log10(ReadNormNum)") +
    theme_bw() +
    theme(text = element_text(family="ArialMT", color = "black", size = 6),
          axis.line = element_line(color = "black"),
          axis.text = element_text(color = "black"),
          axis.ticks = element_line(color = "black"),
          legend.background = element_blank(),
          legend.title = element_text(family="ArialMT", color = "black", size = 6),
          legend.key.size = unit(4, "mm"),
          panel.grid = element_blank())
  
  p1 <- p + theme(legend.position = c(0.8, 0.3))
  ggsave(file.path(oION, sprintf("Nanopore.JS.ht.%s.pdf", samp)), p1, width = 9*2, height = 8, units = "cm")
  
  p2 <- p + theme(legend.position = "none") + geom_abline(slope = 1, intercept = 0)
  ggsave(file.path(oION, sprintf("Nanopore.JS.ht.%s.tiff", samp)), p2, width = 9*2, height = 8, units = "cm", dpi=600)
  
  
  }

write_tsv(all_consistent_gap, "analysis/subgenome/other_data/consistent.SourceData.tsv")

################################################
## Figure S8A

nanopore_bw <- import(file.path(data_dir, "bw", "NarryKim", sprintf("Nanopore.%s.Forward.bw", "VeroInf24h")), "bw")
NGS_bw <- import(file.path("data", "subgenome", "NarryKim_bw", "NarryKim.bw"), "bw")
position_df <- data.frame(Position=1:chrom_size)

signal_NGS <- data.frame(Position=NGS_bw@ranges@start, Signal=NGS_bw$score)
signal_NGS <- left_join(position_df, signal_NGS)
signal_NGS$Source <- "NGS"
signal_Nanopore <- data.frame(Position=nanopore_bw@ranges@start, Signal=nanopore_bw$score)
signal_Nanopore <- left_join(position_df, signal_Nanopore)
signal_Nanopore$Source <- "Nanopore"
merge_signal_df <- rbind(signal_NGS, signal_Nanopore)
merge_signal_df$Signal[is.na(merge_signal_df$Signal)] <- 0

p <- ggplot(merge_signal_df, aes(x=Position, y=log10(Signal+1), color=Source)) +
  geom_line(size=0.1) +
  theme_bw() +
  scale_y_continuous(breaks = log10(c(0+1, 1+1, 10+1, 100+1, 1000+1, 10000+1)), labels = c("0", "1E0", "1E1", "1E2", "1E3", "1E4")) +
  labs(y="Norm. signal (SARS-CoV-2)") +
  theme(text = element_text(family="ArialMT", color = "black", size = 6),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.background = element_blank(),
        legend.title = element_text(family="ArialMT", color = "black", size = 6),
        legend.key.size = unit(4, "mm"),
        panel.grid = element_blank())

ggsave(file.path(file.path(out_dir, "Nanopore", "NarryKim"), sprintf("Nanopore.singal.profile.%s.pdf", "VeroInf24h")), p, width = 7, height = 3.5, units = "cm")

all_consistent_gap <- read_delim("analysis/subgenome/NGS_Nanopore/all_consistent_gap.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
all_consistent_gap <- all_consistent_gap[,1:2]
NGS_JS <- read_delim(sprintf("data/RNA_seq_js/NarryKim/NarryKim_jd.csv", s), "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
NGS_JS <- NGS_JS[, c(4, 5, 8, 16)]
names(NGS_JS) <- c("Start", "End", "Strand", "ReadNum")
NGS_JS <- NGS_JS[NGS_JS$Strand == "+", ]
NGS_JS$Strand <- NULL
names(NGS_JS)[3] <- "NarryKim"
all_consistent_gap <- left_join(all_consistent_gap, NGS_JS)
write_tsv(all_consistent_gap, "analysis/subgenome/other_data/consistent.RNA_seq.SourceData.tsv")
