library(readr)
library(dplyr)
library(ggplot2)
library(rtracklayer)

all_consistent_gap <- read_delim("analysis/subgenome/NGS_Nanopore/all_consistent_gap.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
all_consistent_gap <- all_consistent_gap[,1:2]
for (s in c("6h", "12h", "24h", "48h")){
  NGS_JS <- read_delim(sprintf("data/RNA_seq_js/TP/%s_jd.csv", s), "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  NGS_JS <- NGS_JS[, c(4, 5, 8, 16)]
  names(NGS_JS) <- c("Start", "End", "Strand", "ReadNum")
  NGS_JS <- NGS_JS[NGS_JS$Strand == "+", ]
  NGS_JS$Strand <- NULL
  names(NGS_JS)[3] <- s
  all_consistent_gap <- left_join(all_consistent_gap, NGS_JS)
}
write_tsv(all_consistent_gap, "analysis/subgenome/timepoint/consistent.RNA_seq.SourceData.tsv")

all_signam_df <- data.frame()
for (s in c("6h", "12h", "24h")) {
  nanopore_bw <- import(file.path("data", "subgenome", "other_data", "bw", "TP", sprintf("Nanopore.Vero%s.Forward.bw", s)), "bw")
  position_df <- data.frame(Position=1:chrom_size)
  signal_Nanopore <- data.frame(Position=nanopore_bw@ranges@start, Signal=nanopore_bw$score)
  signal_Nanopore <- left_join(position_df, signal_Nanopore)
  signal_Nanopore$TP <- s
  signal_Nanopore$Signal[is.na(signal_Nanopore$Signal)] <- 0
  all_signam_df <- rbind(all_signam_df, signal_Nanopore)
}

nanopore_bw <- import(file.path("data", "subgenome", "bw", "Nanopore.Forward.bw"), "bw")
position_df <- data.frame(Position=1:chrom_size)
signal_Nanopore <- data.frame(Position=nanopore_bw@ranges@start, Signal=nanopore_bw$score)
signal_Nanopore <- left_join(position_df, signal_Nanopore)
signal_Nanopore$TP <- "48h"
signal_Nanopore$Signal[is.na(signal_Nanopore$Signal)] <- 0
all_signam_df <- rbind(all_signam_df, signal_Nanopore)
all_signam_df$TP <- factor(all_signam_df$TP, levels = c("6h", "12h", "24h", "48h"))

################################################
## Figure 5A

p <- ggplot(all_signam_df, aes(x=Position, y=log10(Signal+1), color=TP)) +
  geom_line(size=0.1, alpha=0.6) +
  theme_bw() +
  scale_y_continuous(breaks = log10(c(0+1, 1+1, 10+1, 100+1, 1000+1, 10000+1)), labels = c("0", "1E0", "1E1", "1E2", "1E3", "1E4")) +
  labs(y="Norm. signal (SARS-CoV-2)") +
  theme(text = element_text(family="ArialMT", color = "black", size = 5),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(family="ArialMT", color = "black", size = 5),
        legend.position = c(0.5, 0.9),
        legend.key.size = unit(2, "mm"),
        panel.grid = element_blank())
ggsave(sprintf("analysis/subgenome/timepoint/Nanopore.singal.profile.all.pdf", s), p, width = 5, height = 3.5, units = "cm")


compute_size_factor <- function(x, bw, window_size = 100){
  return(sum(bw$score[abs(bw@ranges@start-x) < window_size], na.rm = T))
}

#' Junction gap cutoff
min_len_cutoff <- 100

################################################
## Figure 5A

load_NGS_js <- function(f, label, NGS_bw, window_size=100){
  NGS_JS <- read_delim(f, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  NGS_JS <- NGS_JS[, c(4, 5, 8, 15, 16)]
  names(NGS_JS) <- c("Start", "End", "Strand", "EventNum", "ReadNum")
  NGS_JS <- NGS_JS[NGS_JS$Strand == "+", ]
  NGS_JS$StartFactor <- sapply(NGS_JS$Start, compute_size_factor, bw=NGS_bw, window_size=window_size) + 0.01
  NGS_JS$EndFactor <- sapply(NGS_JS$End, compute_size_factor, bw=NGS_bw, window_size=window_size) + 0.01
  NGS_JS$Factor <- sqrt(NGS_JS$StartFactor * NGS_JS$EndFactor) # genometic average
  NGS_JS$EventNormNum <- NGS_JS$EventNum / NGS_JS$Factor
  NGS_JS$EventNormNum <- NGS_JS$EventNormNum / median(NGS_JS$EventNormNum) # rescale
  NGS_JS$ReadNormNum <- NGS_JS$ReadNum / NGS_JS$Factor
  NGS_JS$ReadNormNum <- NGS_JS$ReadNormNum / median(NGS_JS$ReadNormNum)
  NGS_JS$Len <- NGS_JS$End - NGS_JS$Start
  NGS_JS$EventRatio <- NGS_JS$EventNum / sum(NGS_JS$EventNum) # relative ratio
  NGS_JS$ReadRatio <- NGS_JS$ReadNum / sum(NGS_JS$ReadNum)
  NGS_JS$Source <- label
  return(NGS_JS)
}


################################################
## Figure 1G

for (s in c("6h", "12h", "24h", "48h")) {
  NGS_bw <- import(file.path("data", "RNA_TP_bw", sprintf("%s.Forward.bw", s)), "bw")
  f_js <- file.path("data", "RNA_seq_js", "TP", sprintf("%s_jd.csv", s))
  NGS_js <- load_NGS_js(f_js, s, NGS_bw)
  NGS_read_ratio_cutoff <- quantile(NGS_js$ReadRatio, 0.975) 
  NGS_js$ReadPass <- NGS_js$ReadRatio > NGS_read_ratio_cutoff
  NGS_norm_cutoff <- quantile(NGS_js$ReadNormNum[!NGS_js$ReadPass], 0.99)  # Top 1%
  NGS_js$ReadPass <- factor(NGS_js$ReadPass, levels = c(TRUE, FALSE), labels = c("High expr.", "Low expr."))
  NGS_js$ReadPass <- NULL
  enriched_NGS_js <- NGS_js[
    ((NGS_js$ReadRatio > NGS_read_ratio_cutoff) & 
       (NGS_js$ReadNormNum > NGS_norm_cutoff) & 
       (NGS_js$Len > min_len_cutoff)), ]
  genome_size <- 29891
  p <- ggplot() +
    geom_point(data=NGS_js, mapping=aes(x=Start, y=End), color="#92C5DE", alpha=0.3, size=0.1) +
    geom_point(data=enriched_NGS_js, mapping=aes(x=Start, y=End, alpha=log10(ReadNormNum)), size=1, color="red") +
    lims(x=c(0, genome_size), y=c(0, genome_size)) +
    #geom_abline(slope = 1, intercept = 0) +
    labs(x="Start position", y="End position", alpha="log10(ReadNormNum)") +
    theme_bw() +
    theme(text = element_text(family="ArialMT", color = "black", size = 6),
          axis.line = element_line(color = "black"),
          axis.text = element_text(color = "black"),
          axis.ticks = element_line(color = "black"),
          legend.position = "none",
          panel.grid = element_blank())
  ggsave(file.path("analysis/subgenome/timepoint", sprintf("NGS.JS.%s.ht.tiff", s)), p, width = 9, height = 8, units = "cm", dpi=600)
  
  p <- p + geom_abline(slope = 1, intercept = 0)
  ggsave(file.path("analysis/subgenome/timepoint", sprintf("NGS.JS.%s.ht.pdf", s)), p, width = 9, height = 8, units = "cm")
}
