library(readr)
library(dplyr)
library(ggplot2)
setwd("E:/nCoV_VERO")

################################################

## Figure S1B

fl_cutoff <- 45 # 65 - 20, criteria for being full-(sub)genome at both ends
chrom_size <- 29891
NB_start <- 29090
NB_end <- 29870
load_nanopore_len <- function(f_bed, fl=FALSE, filter_by_probe=TRUE, filter_by_overlap_len=NA){
  bed <- read_delim(f_bed, "\t", escape_double = FALSE, col_names = FALSE, col_types = cols(X11 = col_character()), trim_ws = TRUE)
  if(filter_by_probe){
    last_block_len <- sapply(strsplit(bed$X11, ","), function(x){return(as.integer(x[length(x)]))})
    last_block_start <- sapply(strsplit(bed$X12, ","), function(x){return(as.integer(x[length(x)]))})
    last_block_end <- last_block_start + last_block_len
    bed <- bed[(last_block_start<NB_start) & (last_block_end>NB_end),]
  }
  if(fl){
    bed <- bed[(bed$X2 < fl_cutoff) & (bed$X3 > (chrom_size-fl_cutoff)),]
  }
  if(!is.na(filter_by_overlap_len)){
    last_block_len <- sapply(strsplit(bed$X11, ","), function(x){return(as.integer(x[length(x)]))})
    last_block_start <- sapply(strsplit(bed$X12, ","), function(x){return(as.integer(x[length(x)]))})
    last_block_end <- last_block_start + last_block_len
    overlap_df <- data.frame(block_start=last_block_start, block_end=last_block_end, probe_start=NB_start, probe_end=NB_end)
    overlap_df$max_start <- apply(overlap_df[,c("block_start", "probe_start")], 1, max)
    overlap_df$min_end <- apply(overlap_df[,c("block_end", "probe_end")], 1, min)
    overlap_df$overlap_len <- overlap_df$min_end - overlap_df$max_start
    bed <- bed[overlap_df$overlap_len>filter_by_overlap_len,]
  }
  read_len <- sapply(strsplit(bed$X11, ","), function(x){return(sum(as.integer(x)))}) 
  return(read_len)
}

subgenome_group <- read_delim("analysis/subgenome/subgenome/subgenome.group.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
max_subgenome <- subgenome_group %>% group_by(Type) %>% arrange(-GroupReadNumber) %>% slice(1:1)
max_subgenome <- max_subgenome[max_subgenome$Type!="Other",]
js_end <- sapply(strsplit(max_subgenome$LeaderJuncList, "-"), function(x){return(as.integer(x[2]))})
js_start <- sapply(strsplit(max_subgenome$LeaderJuncList, "-"), function(x){return(as.integer(x[1]))})
max_subgenome$SgLen <- chrom_size - (js_end-js_start)

plot_NB <- function(read_len){
  k <- 4.5
  b <- 1
  bw <- 0.01
  filter_len <- 100
  marker_li <- c(1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 9000)
  oe_li <- c(0.1, 0.2, 0.5, 0.7, 0.8, 0.9, 0.95, 0.97, 0.98, 0.99)
  ue_li <- c(0.1, 0.2, 0.5, 0.7)
  
  sgRNA_loc <- data.frame(
    Type=c("S", "NS3", "E", "M", "NS6", "NS7a", "NS7b", "NS8", "N"), 
    Start=c(65, 65, 69, 64, 69, 66, 66, 65, 64),
    End=c(21551, 25380, 26236, 26467, 27040, 27384, 27484, 27883, 28254)
  )
  sgRNA_loc$SgLen <- chrom_size - (sgRNA_loc$End-sgRNA_loc$Start)
  
  
  # lg(M) = -bm + k
  # M: Molecular weight
  # m: Mobility
  
  # m=(k-lg(M))/b
  # m=(k'-log10(l))/b'
  m <- (k - log10(read_len)) / b
  marker_m <- (k - log10(marker_li)) / b
  ref_m <- (k - log10(sgRNA_loc$SgLen)) / b
  
  dense_li <- density(m, n=2000, from=0, to=max(m), bw=bw)
  tmp_df <- data.frame(Type="Best", x=0, m=dense_li$x, expr=dense_li$y)
  backup_df <- tmp_df
  
  for(indx in 1:length(oe_li)){
    oe <- oe_li[indx]
    tmp_oe_df <- backup_df
    max_expr <- max(tmp_oe_df$expr)
    expr_cutoff <- (1-oe) * max_expr
    tmp_oe_df$expr[tmp_oe_df$expr>expr_cutoff] <- expr_cutoff
    tmp_oe_df$expr <- tmp_oe_df$expr / expr_cutoff * max_expr
    tmp_oe_df$Type <- sprintf("Overexposure %d%%", oe*100)
    tmp_oe_df$x <- indx
    tmp_df <- rbind(tmp_df, tmp_oe_df)
  }
  
  for(indx in 1:length(ue_li)){
    ue <- ue_li[indx]
    tmp_ue_df <- backup_df
    tmp_ue_df$expr <- (1-ue) * tmp_ue_df$expr
    tmp_ue_df$Type <- sprintf("Underexposure %d%%", ue*100)
    tmp_ue_df$x <- indx + length(oe_li)
    tmp_df <- rbind(tmp_df, tmp_ue_df)
  }
  
  x_df <- tmp_df %>% group_by(x, Type) %>% summarise()
  tmp_resolution <- tmp_df$m[2] - tmp_df$m[1]
  p <- ggplot() +
    geom_rect(data=tmp_df, mapping = aes(xmin=-0.4+x, xmax=0.4+x, ymin=m-(tmp_resolution/2), ymax=m+(tmp_resolution/2), alpha=expr), fill="black") +
    scale_alpha_continuous(range = c(0, 1)) +
    scale_x_continuous(breaks = x_df$x, labels = x_df$Type) +
    scale_y_reverse(breaks = marker_m, labels = marker_li, sec.axis=sec_axis(~., breaks = ref_m, labels = sgRNA_loc$Type)) +
    theme_bw() +
    theme(text = element_text(family="ArialMT", color = "black", size = 6),
          axis.line = element_line(color = "black"),
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.ticks = element_line(color = "black"),
          legend.position = "none",
          panel.grid = element_blank())
  return(p)
}

rep1_len <- load_nanopore_len(sprintf("data/subgenome/nanopore/data/WIV04.%s.bed", "rep1"))
rep2_len <- load_nanopore_len(sprintf("data/subgenome/nanopore/data/WIV04.%s.bed", "rep2"))
read_len_48h <- c(rep1_len, rep2_len)
p_48h <- plot_NB(read_len_48h)
ggsave("analysis/NB/simulate.48h.pdf", p_48h, width = 9, height = 6, units = "cm")

read_len_12h <- load_nanopore_len("data/subgenome/other_data/nanopore/data/TP/WIV04.Vero12h.bed")
p_12h <- plot_NB(read_len_12h)
ggsave("analysis/NB/simulate.12h.pdf", p_12h, width = 9, height = 6, units = "cm")

read_len_24h <- load_nanopore_len("data/subgenome/other_data/nanopore/data/TP/WIV04.Vero24h.bed")
p_24h <- plot_NB(read_len_24h)
ggsave("analysis/NB/simulate.24h.pdf", p_24h, width = 9, height = 6, units = "cm")
