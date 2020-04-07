# 1. Identify significant junctions

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
if (!dir.exists(out_dir)) dir.create(out_dir)

##### NGS data

#' compute positional signal
compute_size_factor <- function(x, bw, window_size = 100){
  return(sum(bw$score[abs(bw@ranges@start-x) < window_size], na.rm = T))
}

#' Junction gap cutoff
min_len_cutoff <- 100

## Identify NGS JS
js_data_dir <- file.path(data_dir, "js")
NGS_bw <- import(file.path(data_dir, "bw", "NGS.Forward.bw"), "bw")

#' Read and normalize signal
#' f: junctions site detailed information from rnatk_ja
load_NGS_js <- function(f, label, window_size=100){
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

NGS_JS_rep1 <- load_NGS_js(file.path(js_data_dir, "CoV_rep1_jd.csv"), "Rep1")
NGS_JS_rep2 <- load_NGS_js(file.path(js_data_dir, "CoV_rep2_jd.csv"), "Rep2")
NGS_JS_merged <- load_NGS_js(file.path(js_data_dir, "CoV_merged_jd.csv"), "Merged")

all_NGS_JS_df <- rbind(NGS_JS_rep1, NGS_JS_rep2, NGS_JS_merged)
all_NGS_JS_df$Source <- factor(all_NGS_JS_df$Source, levels = c("Rep1", "Rep2", "Merged"))
# Save all junction site
oNGS <- file.path(out_dir, "NGS")
if (!dir.exists(oNGS)) dir.create(oNGS)
write_tsv(all_NGS_JS_df, file.path(oNGS, "cutoff.sourceData.tsv"))

# Cutoff of junction counts normalized to total# of junctions
NGS_read_ratio_cutoff <- quantile(all_NGS_JS_df$ReadRatio, 0.975)  # Top 2.5%

head(all_NGS_JS_df)
summary(all_NGS_JS_df$ReadRatio)
p <- ggplot(all_NGS_JS_df, aes(x=ReadRatio * 1e5, fill=Source)) +
  geom_histogram(bins=25, alpha=0.3) +
  geom_vline(xintercept = NGS_read_ratio_cutoff * 1e5) +
  annotate("text", x= NGS_read_ratio_cutoff * 1e5, y=5000, 
           label=sprintf("Cutoff = %.2e", NGS_read_ratio_cutoff), size=2, color="black", hjust=-0.1) +
  scale_y_log10() +
  facet_wrap(~Source) +
  scale_x_continuous(limits = c(-0.2, 5)) +
  labs(x="Read ratio (x1e-5)", y="#Junctions") +
  theme_bw() +
  theme(text = element_text(family="ArialMT", color = "black", size = 6),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none",
        panel.grid = element_blank())
ggsave(file.path(oNGS, "ReadRatio.cutoff.pdf"), p, width = 12, height = 4, units = "cm")

# Normalize to local background
# Using junctions without enough #reads as background, choose alpha=0.01 to compute cutoff of local background
all_NGS_JS_df$ReadPass <- all_NGS_JS_df$ReadRatio > NGS_read_ratio_cutoff
NGS_norm_cutoff <- quantile(all_NGS_JS_df$ReadNormNum[!all_NGS_JS_df$ReadPass], 0.99)  # Top 1%
all_NGS_JS_df$ReadPass <- factor(all_NGS_JS_df$ReadPass, levels = c(TRUE, FALSE), labels = c("High expr.", "Low expr."))
p <- ggplot(all_NGS_JS_df, aes(x=ReadNormNum, fill=Source)) +
  geom_histogram(bins=50, alpha=0.3, position="identity") +
  geom_vline(xintercept = NGS_norm_cutoff) +
  facet_grid(vars(ReadPass), vars(Source), scales="free_y") +
  scale_x_log10() +
  labs(x="Normalized number", y="#Junctions") +
  theme_bw() +
  theme(text = element_text(family="ArialMT", color = "black", size = 6),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(4, "mm"),
        legend.position = "none",
        panel.grid = element_blank())
ggsave(file.path(oNGS, "ReadNormNum.cutoff.pdf"), p, width = 12, height = 6, units = "cm")

# Save enriched junctions
all_NGS_JS_df$ReadPass <- NULL
enriched_NGS_JS <- all_NGS_JS_df[
  ((all_NGS_JS_df$ReadRatio > NGS_read_ratio_cutoff) & 
   (all_NGS_JS_df$ReadNormNum > NGS_norm_cutoff) & 
   (all_NGS_JS_df$Len > min_len_cutoff)), ]
write_tsv(enriched_NGS_JS, file.path(oNGS, "NGS.enriched.sourceData.tsv"))

# Compare between replicates
enriched_NGS_JS_rep1 <- enriched_NGS_JS[enriched_NGS_JS$Source=="Rep1", ]
enriched_NGS_JS_rep2 <- enriched_NGS_JS[enriched_NGS_JS$Source=="Rep2", ]
enriched_NGS_JS_merged <- enriched_NGS_JS[enriched_NGS_JS$Source=="Merged", ]
enriched_NGS_JS_merged$Source <- NULL

enriched_NGS_JS_merged$Rep1 <- apply(enriched_NGS_JS_merged, 1, function(line, JS_df) {
  res <- any((JS_df$Start == as.integer(line["Start"])) & (JS_df$End == as.integer(line["End"])))
  return(res)
}, JS_df = enriched_NGS_JS_rep1)

enriched_NGS_JS_merged$Rep2 <- apply(enriched_NGS_JS_merged, 1, function(line, JS_df){
  res <- any((JS_df$Start == as.integer(line["Start"])) & (JS_df$End == as.integer(line["End"])))
  return(res)
}, JS_df = enriched_NGS_JS_rep2)

enriched_NGS_JS_merged$Consistent <- enriched_NGS_JS_merged$Rep1 & enriched_NGS_JS_merged$Rep2
enriched_NGS_JS_merged <- enriched_NGS_JS_merged[order(enriched_NGS_JS_merged$ReadNormNum, decreasing = T), ]
write_tsv(enriched_NGS_JS_merged, file.path(oNGS, "NGS.enriched.tsv"))
# save reproducible junctions
enriched_NGS_JS_merged <- enriched_NGS_JS_merged[enriched_NGS_JS_merged$Consistent, ]
write_tsv(enriched_NGS_JS_merged, file.path(oNGS, "NGS.enriched.consistent.tsv"))

#' Merge neighbouring junctions within window into group
merge_gap <- function(gap_df, resolution = 5) {
  gap_df <- gap_df[, c("Start", "End", "ReadNum")]
  # check if two junctions can be merged
  is_one_group <- function(gap1, gap2, resolution = 5) {
    diff_gap <- gap1 - gap2
    return(all(abs(diff_gap) <= resolution))
  }
  indx1 <- c()
  indx2 <- c()
  for (i in 1:nrow(gap_df)) {
    for (j in 1:i) {
      can_merge <- is_one_group(gap_df[i, c("Start", "End")], gap_df[j, c("Start", "End")], resolution)
      if (can_merge) {
        indx1 <- c(indx1, i)
        indx2 <- c(indx2, j)
      }
    }
  }
  edge_df <- data.frame(from = indx1, to = indx2)
  edge_df <- edge_df[order(edge_df$from, edge_df$to), ]
  gap_graph <- graph_from_edgelist(as.matrix(edge_df), directed = FALSE)
  graph_clustar <- components(gap_graph, "strong")
  
  merge_gap <- function(indx_li, gap_df) {
    selected_df <- gap_df[indx_li, ]
    selected_df <- selected_df[order(selected_df$ReadNum, decreasing = TRUE), ] 
    res <- data.frame(
      LeaderStart = selected_df$Start[1], # Core junction with maximum reads in one group
      LeaderEnd = selected_df$End[1],
      GroupNum = nrow(selected_df),
      GroupReadNum = sum(selected_df$ReadNum)
    )
    return(res)
  }
  merged_group <- do.call(rbind, lapply(groups(graph_clustar), merge_gap, gap_df = gap_df))
  merged_group <- merged_group[order(merged_group$GroupReadNum, decreasing = TRUE), ]
  return(merged_group)
}

NGS_JS_group <-  merge_gap(enriched_NGS_JS_merged)
# save junction group 
write_tsv(NGS_JS_group, file.path(oNGS, "consistent_gap.group.tsv"))

# Show all NGS junctions and highight enriched ones
genome_size <- 29891
p <- ggplot() +
  geom_point(data=NGS_JS_merged, mapping=aes(x=Start, y=End), color="grey80", alpha=0.3, size=0.1) +
  geom_point(data=enriched_NGS_JS_merged, mapping=aes(x=Start, y=End, alpha=log10(ReadNormNum)), size=1, color="red") +
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
ggsave(file.path(oNGS, "NGS.JS.ht.tiff"), p, width = 9, height = 8, units = "cm", dpi=600)

p <- p + geom_abline(slope = 1, intercept = 0)
ggsave(file.path(oNGS, "NGS.JS.ht.pdf"), p, width = 9, height = 8, units = "cm")

# Show comparison between replicates
junc_df <- enriched_NGS_JS_merged[, c("Start", "End")]
junc_df$Enriched <- TRUE
rep1_df <- NGS_JS_rep1[, c("Start", "End", "ReadNum")]
rep2_df <- NGS_JS_rep2[, c("Start", "End", "ReadNum")]
names(rep1_df)[3] <- "Rep1"
names(rep2_df)[3] <- "Rep2"
jun_all_df <- full_join(rep1_df, rep2_df)
head(jun_all_df)
jun_all_df[is.na(jun_all_df)] <- 0

jun_all_df <- left_join(jun_all_df, junc_df)
jun_all_df$Enriched[is.na(jun_all_df$Enriched)] <- FALSE
p <- ggplot() +
  geom_point(data=jun_all_df, mapping=aes(x=Rep1+1, y=Rep2+1), color="grey70", alpha=0.3, size=0.4) +
  geom_point(data=jun_all_df[jun_all_df$Enriched, ], mapping=aes(x=Rep1+1, y=Rep2+1), color="red", size=0.5) +
  labs(x="#Reads+1 in Rep1", y="#Reads+1 in Rep2") +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  theme(text = element_text(family="ArialMT", color = "black", size = 6),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        panel.grid = element_blank(),
        panel.border = element_blank()
        )
ggsave(file.path(oNGS, "NGS.JS.cor.pdf"), p, width = 9, height = 8, units = "cm")
ggsave(file.path(oNGS, "NGS.JS.cor.tiff"), p, width = 9, height = 8, units = "cm", dpi=600)



###### 2 Nanopore data

fl_cutoff <- 45 # 65 - 20, criteria for being full-(sub)genome at both ends
chrom_size <- 29891
nanopore_bw <- import(file.path(data_dir, "bw", "Nanopore.Forward.bw"), "bw")
nanopore_bed_rep1 <- read_delim(file.path(data_dir, "nanopore", "data", "WIV04.rep1.bed"), 
                                "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
nanopore_bed_rep2 <- read_delim(file.path(data_dir, "nanopore", "data", "WIV04.rep2.bed"), 
                                "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
# keep only full-length reads for downstream analysis
nanopore_rep1_fl <- nanopore_bed_rep1$X4[(nanopore_bed_rep1$X2 < fl_cutoff) & (nanopore_bed_rep1$X3 > (chrom_size-fl_cutoff))]
nanopore_rep2_fl <- nanopore_bed_rep2$X4[(nanopore_bed_rep2$X2 < fl_cutoff) & (nanopore_bed_rep2$X3 > (chrom_size-fl_cutoff))]

#' Read Nanopore junctions and normalize the signals
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

nanopore_JS_rep1 <- load_nanopore_info(file.path(data_dir, "nanopore", "data", "WIV04.rep1.gap.bed"), nanopore_rep1_fl, "Rep1")
nanopore_JS_rep2 <- load_nanopore_info(file.path(data_dir, "nanopore", "data", "WIV04.rep2.gap.bed"), nanopore_rep2_fl, "Rep2")
# Do same analysis for merged data
nanopore_JS_merged <- rbind(nanopore_JS_rep1, nanopore_JS_rep2)
nanopore_JS_merged <- nanopore_JS_merged %>% group_by(FL, Start, End) %>% summarise(ReadNum=sum(ReadNum))
nanopore_JS_merged$StartFactor <- sapply(nanopore_JS_merged$Start, compute_size_factor, bw=nanopore_bw, window_size=100) + 0.01
nanopore_JS_merged$EndFactor <- sapply(nanopore_JS_merged$End, compute_size_factor, bw=nanopore_bw, window_size=100) + 0.01
nanopore_JS_merged$Factor <- sqrt(nanopore_JS_merged$StartFactor * nanopore_JS_merged$EndFactor)
nanopore_JS_merged$ReadNormNum <- nanopore_JS_merged$ReadNum / nanopore_JS_merged$Factor
nanopore_JS_merged$ReadNormNum <- nanopore_JS_merged$ReadNormNum / median(nanopore_JS_merged$ReadNormNum)
nanopore_JS_merged$Len <- nanopore_JS_merged$End - nanopore_JS_merged$Start
nanopore_JS_merged$ReadRatio <- nanopore_JS_merged$ReadNum / sum(nanopore_JS_merged$ReadNum)
nanopore_JS_merged <- nanopore_JS_merged %>% group_by(FL) %>% mutate(ReadRatioGroup=ReadNum/sum(ReadNum))
nanopore_JS_merged$Source <- "Merged"

all_nanopore_JS_df <- rbind(nanopore_JS_rep1, nanopore_JS_rep2, nanopore_JS_merged)
all_nanopore_JS_df$Source <- factor(all_nanopore_JS_df$Source, levels = c("Rep1", "Rep2", "Merged"))
oION <- file.path(out_dir, "Nanopore")
if (!dir.exists(oION)) dir.create(oION)
write_tsv(all_nanopore_JS_df, file.path(oION, "cutoff.sourceData.tsv"))

# Require more than 2 Nanopore reads for valid junctions
nanopore_read_cutoff <- 2
all_nanopore_JS_df$ReadPass <- all_nanopore_JS_df$ReadNum > nanopore_read_cutoff

p <- ggplot(all_nanopore_JS_df, aes(x=ReadNum, fill=ReadPass)) +
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
ggsave(file.path(oION, "ReadNum.cutoff.pdf"), p, width = 12, height = 6, units = "cm")

# Compute Nanopore local background cutoff
nanopore_norm_cutoff <- quantile(all_nanopore_JS_df$ReadNormNum[!all_nanopore_JS_df$ReadPass], 0.99)  # Top 1%
all_nanopore_JS_df$ReadPass <- factor(all_nanopore_JS_df$ReadPass, levels = c(TRUE, FALSE), labels = c("High expr.", "Low expr."))
p <- ggplot(all_nanopore_JS_df, aes(x=ReadNormNum, fill=Source)) +
  geom_histogram(bins=50, alpha=0.3, position="identity") +
  geom_vline(xintercept = nanopore_norm_cutoff, size=0.3, lty="dotted") +
  facet_grid(vars(ReadPass, FL), vars(Source), scales="free_y") +
  scale_x_log10() +
  labs(x="Normalized number", y="#Junctions") +
  theme_bw() +
  theme(text = element_text(family="ArialMT", color = "black", size = 6),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none",
        panel.grid = element_blank())
ggsave(file.path(oION, "ReadNormNum.cutoff.pdf"), p, width = 12, height = 8, units = "cm")

# Identify Nanopore junctions passing local background and read# cutoff 
all_nanopore_JS_df$ReadPass <- NULL
enriched_nanopore_JS <- all_nanopore_JS_df[
  ((all_nanopore_JS_df$ReadNum>nanopore_read_cutoff) & 
   (all_nanopore_JS_df$ReadNormNum>nanopore_norm_cutoff) & 
   (all_nanopore_JS_df$Len>min_len_cutoff)), ]
write_tsv(enriched_nanopore_JS, file.path(oION, "Nanopore.enriched.sourceData.tsv"))

enriched_nanopore_JS_rep1 <- enriched_nanopore_JS[enriched_nanopore_JS$Source=="Rep1", ]
enriched_nanopore_JS_rep2 <- enriched_nanopore_JS[enriched_nanopore_JS$Source=="Rep2", ]
enriched_nanopore_JS_merged <- enriched_nanopore_JS[enriched_nanopore_JS$Source=="Merged", ]
enriched_nanopore_JS_merged$Source <- NULL

enriched_nanopore_JS_merged$Rep1 <- apply(enriched_nanopore_JS_merged, 1, function(line, JS_df){
  res <- any((JS_df$Start == as.integer(line["Start"])) & (JS_df$End == as.integer(line["End"]) & (JS_df$FL == line["FL"])))
  return(res)
}, JS_df=enriched_nanopore_JS_rep1)

enriched_nanopore_JS_merged$Rep2 <- apply(enriched_nanopore_JS_merged, 1, function(line, JS_df){
  res <- any((JS_df$Start == as.integer(line["Start"])) & (JS_df$End == as.integer(line["End"]) & (JS_df$FL == line["FL"])))
  return(res)
}, JS_df=enriched_nanopore_JS_rep2)

enriched_nanopore_JS_merged$Consistent <- enriched_nanopore_JS_merged$Rep1 & enriched_nanopore_JS_merged$Rep2
enriched_nanopore_JS_merged <- enriched_nanopore_JS_merged[
  order(enriched_nanopore_JS_merged$FL, enriched_nanopore_JS_merged$ReadNormNum, decreasing = T), ]
write_tsv(enriched_nanopore_JS_merged, file.path(oION, "Nanopore.enriched.tsv"))

enriched_nanopore_JS_merged <- enriched_nanopore_JS_merged[enriched_nanopore_JS_merged$Consistent,]
write_tsv(enriched_nanopore_JS_merged, file.path(oION, "Nanopore.enriched.consistent.tsv"))

# Show reproduicity between replicates
p <- ggplot() +
  geom_point(data=nanopore_JS_merged, mapping=aes(x=Start, y=End), color="grey80", size=0.1, alpha=0.3) +
  geom_point(data=enriched_nanopore_JS_merged, mapping=aes(x=Start, y=End, alpha=log10(ReadNormNum)), size=1, color="red") +
  lims(x=c(0, genome_size), y=c(0, genome_size)) +
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
ggsave(file.path(oION, "Nanopore.JS.ht.pdf"), p1, width = 9*2, height = 8, units = "cm")

p2 <- p + theme(legend.position = "none") + geom_abline(slope = 1, intercept = 0)
ggsave(file.path(oION, "Nanopore.JS.ht.tiff"), p2, width = 9*2, height = 8, units = "cm", dpi=600)

# Compare Nanopore replicates
junc_df <- enriched_nanopore_JS_merged[, c("FL", "Start", "End")]
junc_df$Enriched <- TRUE
rep1_df <- nanopore_JS_rep1[, c("FL", "Start", "End", "ReadNum")]
rep2_df <- nanopore_JS_rep2[, c("FL", "Start", "End", "ReadNum")]
names(rep1_df)[4] <- "Rep1"
names(rep2_df)[4] <- "Rep2"
jun_all_df <- full_join(rep1_df, rep2_df)
jun_all_df[is.na(jun_all_df)] <- 0
jun_all_df <- left_join(jun_all_df, junc_df)
jun_all_df$Enriched[is.na(jun_all_df$Enriched)] <- FALSE
p <- ggplot() +
  geom_point(data=jun_all_df, mapping=aes(x=Rep1+1, y=Rep2+1), color="grey70", alpha=0.1, size=0.4) +
  geom_point(data=jun_all_df[jun_all_df$Enriched, ], mapping=aes(x=Rep1+1, y=Rep2+1), color="red", size=0.5) +
  labs(x="#Reads+1 in Rep1", y="#Reads+1 in Rep2") +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~FL) +
  theme_bw() +
  theme(text = element_text(family="ArialMT", color = "black", size = 6),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        panel.grid = element_blank())
ggsave(file.path(oION, "Nanopore.JS.cor.pdf"), p, width = 7.5*2, height = 7.5, units = "cm")
ggsave(file.path(oION, "Nanopore.JS.cor.tiff"), p, width = 7.5*2, height = 7.5, units = "cm", dpi=600)


###### Compare NGS and Nanopore data
enriched_nanopore_JS_merged_FL <- enriched_nanopore_JS_merged[enriched_nanopore_JS_merged$FL=="FL", ]
enriched_nanopore_JS_merged_NFL <- enriched_nanopore_JS_merged[enriched_nanopore_JS_merged$FL=="non-FL", ]
junc_li <- list(
  "RNA-seq" = paste(enriched_NGS_JS_merged$Start, enriched_NGS_JS_merged$End),
  "Nanopore (FL)" = paste(enriched_nanopore_JS_merged_FL$Start, enriched_nanopore_JS_merged_FL$End),
  "Nanopore (non-FL)" = paste(enriched_nanopore_JS_merged_NFL$Start, enriched_nanopore_JS_merged_NFL$End)
)

oNGSION <- file.path(out_dir, "NGS_Nanopore")
if (!dir.exists(oNGSION)) dir.create(oNGSION)
venn.diagram(
  junc_li,
  file.path(out_dir, "NGS_Nanopore", "NGS_Nanopore.enrich.overlap.venn.tiff"),
  fill = brewer.pal(3, "Set1"),
  height = 450*4, width = 450*4, resolution = 150*4,
  cat.cex = 0.55
)

# save data for the venn plot
NGS_junc_df <- enriched_NGS_JS_merged[, c("Start", "End")]
NGS_junc_df$NGS <- TRUE
nanopore_junc_FL_df <- enriched_nanopore_JS_merged_FL[, c("Start", "End")]
nanopore_junc_FL_df$Nanopore_FL <- TRUE
nanopore_junc_NFL_df <- enriched_nanopore_JS_merged_NFL[, c("Start", "End")]
nanopore_junc_NFL_df$Nanopore_nFL <- TRUE
junc_df <- full_join(NGS_junc_df, nanopore_junc_FL_df)
junc_df <- full_join(junc_df, nanopore_junc_NFL_df)
junc_df[is.na(junc_df)] <- FALSE
write_tsv(junc_df, file.path(out_dir, "NGS_Nanopore", "NGS_Nanopore.enrich.overlap.sourceData.tsv"))


junc_df$Tag <- paste(junc_df$NGS, junc_df$Nanopore_FL)
junc_df$Tag <- factor(junc_df$Tag, 
                      levels = c("TRUE TRUE", "TRUE FALSE", "FALSE TRUE"), 
                      labels = c("Both enriched", "NGS enriched", "Nanopore enriched"))

p <- ggplot(junc_df, aes(x=Start, y=End, color=Tag, size=Tag)) +
  geom_point(alpha = 0.4) +
  lims(x=c(0, genome_size), y=c(0, genome_size)) +
  geom_abline(slope = 1, intercept = 0) +
  scale_size_manual(values = c("Both enriched"=1.2, "NGS enriched"=0.3, "Nanopore enriched"=0.3)) +
  labs(x="Start position", y="End position") +
  theme_bw() +
  theme(text = element_text(family="ArialMT", color = "black", size = 6),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(4, "mm"),
        legend.position = c(0.8, 0.3),
        panel.grid = element_blank())
ggsave(file.path(oNGSION, "NGS_Nanopore.enrich.pdf"), p, width = 8.5, height = 8, units = "cm")

# Identify junctions supported by both NGS and Nanopore including non-enriched ones in NGS/Nanopre
all_NGS_js <- inner_join(NGS_JS_rep1[, c("Start", "End", "ReadNum")], 
                         NGS_JS_rep2[, c("Start", "End", "ReadNum")], by=c("Start", "End"))
names(all_NGS_js)[3:4] <- c("NGS_rep1_Cnt", "NGS_rep2_Cnt")
all_Nanopore_js <- inner_join(nanopore_JS_rep1[nanopore_JS_rep1$FL == "FL", c("Start", "End", "ReadNum")], 
                              nanopore_JS_rep2[nanopore_JS_rep2$FL == "FL", c("Start", "End", "ReadNum")], by=c("Start", "End"))
names(all_Nanopore_js)[3:4] <- c("Nanopore_rep1_Cnt", "Nanopore_rep2_Cnt")
all_consistent_js <- inner_join(all_NGS_js, all_Nanopore_js)
all_consistent_js <- all_consistent_js[!duplicated(all_consistent_js), ]
write_tsv(all_consistent_js, file.path(oNGSION, "all_consistent_gap.tsv"))
nrow(all_consistent_js)

p <- ggplot() +
  geom_point(data=all_consistent_js, mapping=aes(x=Start, y=End), size=0.3, color="black") +
  lims(x=c(0, genome_size), y=c(0, genome_size)) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x="Start position", y="End position") +
  theme_bw() +
  theme(text = element_text(family="ArialMT", color = "black", size = 6),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.background = element_blank(),
        legend.title = element_text(family="ArialMT", color = "black", size = 6),
        legend.key.size = unit(4, "mm"),
        legend.position = c(0.8, 0.3),
        panel.grid = element_blank())
ggsave(file.path(oNGSION, "all_consistent_gap.pdf"), p, width = 9, height = 8, units = "cm")




