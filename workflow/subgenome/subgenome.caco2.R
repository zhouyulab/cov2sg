library(readr)
library(dplyr)
library(ggplot2)


all_consistent_gap <- read_delim("analysis/subgenome/NGS_Nanopore/all_consistent_gap.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

fl_cutoff <- 45 # 65 - 20, criteria for being full-(sub)genome at both ends
chrom_size <- 29891

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

compute_size_factor <- function(x, bw, window_size = 100){
  return(sum(bw$score[abs(bw@ranges@start-x) < window_size], na.rm = T))
}

load_nanopore_info <- function(f, fl_li, label, nanopore_bw, window_size=100){
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


tps <- c("12h", "24h")
rep = "rep1"

total_reads <- c()
fl_reads <- c()
js_df <- data.frame()

################################################
## Figure S3C

for(tp in tps){
  caco2_sample <- sprintf("%s.%s", tp, rep)
  js_info <- load_Nanopore_JS(sprintf("data/Caco2/nanopore/data/Caco2.%s.bed", caco2_sample), sprintf("data/Caco2/nanopore/data/Caco2.%s.gap.bed", caco2_sample))
  total_reads <- c(total_reads, js_info$total_reads)
  fl_reads <- c(fl_reads, js_info$fl_reads)
  tmp_js_df <- js_info$gap
  tmp_js_df$TP <- tp
  js_df <- rbind(js_df, as.data.frame(tmp_js_df))
  tmp_js_df <- tmp_js_df[,c("Start", "End", "ReadNum")]
  names(tmp_js_df)[3] <- tp
  all_consistent_gap <- left_join(all_consistent_gap, tmp_js_df)
  nanopore_bw <- import(file.path("data", "Caco2", "Nanopore_bw", sprintf("Caco2.%s.rep1.Forward.bw", tp)), "bw")
  nanopore_bed <- read_delim(sprintf("data/Caco2/nanopore/data/Caco2.%s.rep1.bed", tp), "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  nanopore_fl <- nanopore_bed$X4[(nanopore_bed$X2 < fl_cutoff) & (nanopore_bed$X3 > (chrom_size-fl_cutoff))]
  nanopore_JS <- load_nanopore_info(sprintf("data/Caco2/nanopore/data/Caco2.%s.rep1.gap.bed", tp), nanopore_fl, tp, nanopore_bw, window_size=100)
  nanopore_read_cutoff <- 2
  nanopore_JS$ReadPass <- nanopore_JS$ReadNum > nanopore_read_cutoff
  nanopore_norm_cutoff <- quantile(nanopore_JS$ReadNormNum[!nanopore_JS$ReadPass], 0.99)  # Top 1%
  nanopore_JS$Enrich <- nanopore_JS$ReadPass & nanopore_JS$ReadNormNum > nanopore_norm_cutoff
  enriched_nanopore_JS <- nanopore_JS[nanopore_JS$Enrich,]
  p <- ggplot() +
    geom_point(data=nanopore_JS, mapping=aes(x=Start, y=End), color="#92C5DE", size=0.1, alpha=0.3) +
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
  ggsave(file.path("analysis/subgenome/caco2", sprintf("Nanopore.JS.ht.%s.pdf", tp)), p1, width = 8*2, height = 8, units = "cm")
  
  p2 <- p + theme(legend.position = "none") + geom_abline(slope = 1, intercept = 0)
  ggsave(file.path("analysis/subgenome/caco2", sprintf("Nanopore.JS.ht.%s.tiff", tp)), p2, width = 8*2, height = 8, units = "cm", dpi=600)
  
}
write_tsv(all_consistent_gap, "analysis/subgenome/caco2/consistent.SourceData.tsv")

js_df <- js_df[,c("TP", "Start", "End", "ReadNum", "ReadRatio")]
write_tsv(js_df, "analysis/subgenome/caco2/JuncNum.SourceData.tsv")

read_df <- data.frame(TP=tps, TotalReads=total_reads, FullLengthReads=fl_reads)
read_df$FullLengthRatio <- read_df$FullLengthReads / read_df$TotalReads
write_tsv(read_df, "analysis/subgenome/caco2/ReadNum.tsv")
read_df$TP <- factor(read_df$TP, levels = read_df$TP)



all_consistent_gap <- read_delim("analysis/subgenome/NGS_Nanopore/all_consistent_gap.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
all_consistent_gap <- all_consistent_gap[,1:2]
for (s in c("12h", "24h")){
  NGS_JS <- read_delim(sprintf("data/caco2_RNA_seq/js/%s_jd.csv", s), "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  NGS_JS <- NGS_JS[, c(4, 5, 8, 16)]
  names(NGS_JS) <- c("Start", "End", "Strand", "ReadNum")
  NGS_JS <- NGS_JS[NGS_JS$Strand == "+", ]
  NGS_JS$Strand <- NULL
  names(NGS_JS)[3] <- s
  all_consistent_gap <- left_join(all_consistent_gap, NGS_JS)
}
write_tsv(all_consistent_gap, "analysis/subgenome/caco2/consistent.RNA_seq.SourceData.tsv")

################################################
## Figure S8B

all_signam_df <- data.frame()
for (s in c("12h", "24h")) {
  nanopore_bw <- import(file.path("data", "Caco2", "Nanopore_bw", sprintf("Caco2.%s.rep1.Forward.bw", s)), "bw")
  position_df <- data.frame(Position=1:chrom_size)
  signal_Nanopore <- data.frame(Position=nanopore_bw@ranges@start, Signal=nanopore_bw$score)
  signal_Nanopore <- left_join(position_df, signal_Nanopore)
  signal_Nanopore$TP <- s
  signal_Nanopore$Signal[is.na(signal_Nanopore$Signal)] <- 0
  all_signam_df <- rbind(all_signam_df, signal_Nanopore)
}
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
ggsave(sprintf("analysis/subgenome/caco2/Nanopore.singal.profile.all.pdf", s), p, width = 5, height = 3.5, units = "cm")



################################################
## Figure 1I (Caco-2)

JS_cnt_df <- read_delim("analysis/subgenome/caco2/consistent.SourceData.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
JS_cnt_df <- JS_cnt_df[,c("Start", "End", "12h", "24h")]

gene_ref <- read_delim(file.path("data", "subgenome", "genome/WIV04.bed"),
                       "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
gene_ref_df <- data.frame(Name=gene_ref$X4, Start=gene_ref$X2, End=gene_ref$X3)
gene_ref_df <- gene_ref_df[order(gene_ref_df$Start), ]
gene_ref_df$Name <- as.character(gene_ref_df$Name)

leader_js_df <- JS_cnt_df[JS_cnt_df$Start<100,]

leader_js_df$SubgenomeType <- sapply(leader_js_df$End, function (end_site) {
  # locate the first downstream complete gene of the first junction
  idx <- which(end_site < gene_ref_df$Start)
  if (length(idx) == 0) {
    sgtype <- "Other"
  } else {
    sgtype <- gene_ref_df$Name[min(idx)] 
  }
  return(sgtype)}
)

leader_js_cnt_df <- leader_js_df[,c("SubgenomeType", "12h", "24h")]
melt_leader_js_cnt_df <- melt(leader_js_cnt_df, id.vars = "SubgenomeType", variable.name = "Time", value.name = "ReadNum")
melt_leader_js_cnt_df$ReadNum <- as.integer(melt_leader_js_cnt_df$ReadNum)
library_size <- melt_leader_js_cnt_df %>% 
  group_by(Time) %>% 
  summarise(TotalRead=sum(ReadNum, na.rm = TRUE))
subgenome_type_df <- melt_leader_js_cnt_df %>% 
  group_by(SubgenomeType, Time) %>% 
  summarise(ReadNum=sum(ReadNum, na.rm = TRUE))
subgenome_type_df$SubgenomeType <- factor(subgenome_type_df$SubgenomeType, levels = c("Other", gene_ref_df$Name[gene_ref_df$Name!="orf1ab"]))
subgenome_type_df$Time <- factor(subgenome_type_df$Time, levels = c("12h", "24h"))
subgenome_type_df <- left_join(subgenome_type_df, library_size)
subgenome_type_df$ReadRatio <- subgenome_type_df$ReadNum / subgenome_type_df$TotalRead

SubgenomeType_col <- rev(c(brewer.pal(length(levels(subgenome_type_df$SubgenomeType))-1, "RdBu"), "grey70"))
names(SubgenomeType_col) <- levels(subgenome_type_df$SubgenomeType)

p <- ggplot(subgenome_type_df, aes(x=Time, y=ReadRatio, fill=SubgenomeType)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = SubgenomeType_col) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0%", "50%", "100%"), expand = c(0, 0)) +
  labs(y="%sgRNA") +
  theme_bw() +
  theme(text = element_text(family="ArialMT", color = "black", size = 6),
        axis.title.x = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.background = element_blank(),
        legend.title = element_text(family="ArialMT", color = "black", size = 6),
        legend.key.size = unit(3, "mm"),
        panel.grid = element_blank())
ggsave(file.path("analysis/subgenome/caco2", "subgenome.type.cnt.tp.pdf"), p, width = 4.5, height = 4, units = "cm")


NGS_JS_cnt_df <- read_delim("analysis/subgenome/caco2/consistent.RNA_seq.SourceData.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
NGS_JS_cnt_df <- NGS_JS_cnt_df[,c("Start", "End", "12h", "24h")]

NGS_leader_js_df <- NGS_JS_cnt_df[NGS_JS_cnt_df$Start<100,]

NGS_leader_js_df$SubgenomeType <- sapply(NGS_leader_js_df$End, function (end_site) {
  # locate the first downstream complete gene of the first junction
  idx <- which(end_site < gene_ref_df$Start)
  if (length(idx) == 0) {
    sgtype <- "Other"
  } else {
    sgtype <- gene_ref_df$Name[min(idx)] 
  }
  return(sgtype)}
)

NGS_leader_js_cnt_df <- NGS_leader_js_df[,c("SubgenomeType", "12h", "24h")]
NGS_melt_leader_js_cnt_df <- melt(NGS_leader_js_cnt_df, id.vars = "SubgenomeType", variable.name = "Time", value.name = "ReadNum")
NGS_melt_leader_js_cnt_df$ReadNum <- as.integer(NGS_melt_leader_js_cnt_df$ReadNum)
NGS_library_size <- NGS_melt_leader_js_cnt_df %>% 
  group_by(Time) %>% 
  summarise(TotalRead=sum(ReadNum, na.rm = TRUE))
NGS_subgenome_type_df <- NGS_melt_leader_js_cnt_df %>% 
  group_by(SubgenomeType, Time) %>% 
  summarise(ReadNum=sum(ReadNum, na.rm = TRUE))
NGS_subgenome_type_df$SubgenomeType <- factor(NGS_subgenome_type_df$SubgenomeType, levels = c("Other", gene_ref_df$Name[gene_ref_df$Name!="orf1ab"]))
NGS_subgenome_type_df$Time <- factor(NGS_subgenome_type_df$Time, levels = c("6h", "12h", "24h", "48h"))
NGS_subgenome_type_df <- left_join(NGS_subgenome_type_df, NGS_library_size)
NGS_subgenome_type_df$ReadRatio <- NGS_subgenome_type_df$ReadNum / NGS_subgenome_type_df$TotalRead

NGS_SubgenomeType_col <- rev(c(brewer.pal(length(levels(NGS_subgenome_type_df$SubgenomeType))-1, "RdBu"), "grey70"))
names(NGS_SubgenomeType_col) <- levels(NGS_subgenome_type_df$SubgenomeType)

p <- ggplot(NGS_subgenome_type_df, aes(x=Time, y=ReadRatio, fill=SubgenomeType)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = NGS_SubgenomeType_col) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0%", "50%", "100%"), expand = c(0, 0)) +
  labs(y="%sgRNA") +
  theme_bw() +
  theme(text = element_text(family="ArialMT", color = "black", size = 6),
        axis.title.x = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.background = element_blank(),
        legend.title = element_text(family="ArialMT", color = "black", size = 6),
        legend.key.size = unit(3, "mm"),
        panel.grid = element_blank())
p
ggsave(file.path("analysis/subgenome/caco2", "subgenome.type.cnt.tp.NGS.pdf"), p, width = 4.5, height = 4, units = "cm")


compute_size_factor <- function(x, bw, window_size = 100){
  return(sum(bw$score[abs(bw@ranges@start-x) < window_size], na.rm = T))
}

#' Junction gap cutoff
min_len_cutoff <- 100



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
## Figure 1H
for (s in c("12h", "24h")) {
  NGS_bw <- import(file.path("data", "caco2_RNA_seq", "bw", sprintf("%s.Forward.bw", s)), "bw")
  f_js <- file.path("data", "caco2_RNA_seq", "js", sprintf("%s_jd.csv", s))
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
  ggsave(file.path("analysis/subgenome/caco2", sprintf("NGS.JS.%s.ht.tiff", s)), p, width = 9, height = 8, units = "cm", dpi=600)
  
  p <- p + geom_abline(slope = 1, intercept = 0)
  ggsave(file.path("analysis/subgenome/caco2", sprintf("NGS.JS.%s.ht.pdf", s)), p, width = 9, height = 8, units = "cm")
  
}
