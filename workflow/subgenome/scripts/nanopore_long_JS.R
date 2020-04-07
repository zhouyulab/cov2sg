library(readr)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(reshape2)
library(rtracklayer)
library(trackViewer)


data_dir <- file.path("data", "subgenome")
out_dir <- file.path("analysis", "subgenome")
if (!dir.exists(out_dir)) dir.create(out_dir)

chrom_size <- 29891

fl_cutoff <- 45 # 65 - 20, criteria for being full-(sub)genome at both ends

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

# keep only full-length reads for downstream analysis
nanopore_rep1_fl <- nanopore_bed_rep1$X4[(nanopore_bed_rep1$X2 < fl_cutoff) & (nanopore_bed_rep1$X3 > (chrom_size-fl_cutoff))]
nanopore_rep2_fl <- nanopore_bed_rep2$X4[(nanopore_bed_rep2$X2 < fl_cutoff) & (nanopore_bed_rep2$X3 > (chrom_size-fl_cutoff))]
nanopore_multiple_rep1 <- nanopore_bed_rep1 %>% group_by(X4) %>% filter(n()>1) %>% summarise(MapNum=n())
nanopore_multiple_rep2 <- nanopore_bed_rep2 %>% group_by(X4) %>% filter(n()>1) %>% summarise(MapNum=n())
nanopore_rep1_fl <- setdiff(nanopore_rep1_fl, nanopore_multiple_rep1$X4)
nanopore_rep2_fl <- setdiff(nanopore_rep2_fl, nanopore_multiple_rep2$X4)
nanopore_bed_rep1 <- nanopore_bed_rep1[!nanopore_bed_rep1$X4 %in% nanopore_multiple_rep1$X4,]
nanopore_bed_rep2 <- nanopore_bed_rep2[!nanopore_bed_rep2$X4 %in% nanopore_multiple_rep2$X4,]

full_genome <- rbind(
  nanopore_bed_rep1[!nanopore_bed_rep1$X4 %in% nanopore_multiple_rep1$X4,], 
  nanopore_bed_rep2[!nanopore_bed_rep2$X4 %in% nanopore_multiple_rep2$X4,]
  )
full_genome <- full_genome[
  (full_genome$X10==1) &
  (full_genome$X6=="+") &
  (full_genome$X2<fl_cutoff) &
  (full_genome$X3>(chrom_size-fl_cutoff))
  ,]
nrow(full_genome)

load_nanopore_info <- function(f, fl_li, multiple_li, label){
  gap_df <- read_delim(f, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  names(gap_df) <- c("Chrom", "Start", "End", "Name", "Score", "Strand")
  gap_df$FL <- gap_df$Name %in% fl_li
  gap_df <- gap_df[gap_df$Strand == "+", ]
  gap_df <- gap_df[!gap_df$Name %in% multiple_li,]
  gap_df$Chrom <- NULL
  gap_df$Score <- NULL
  gap_df$Strand <- NULL
  gap_df$Source <- label
  return(gap_df)
}

nanopore_JS_rep1 <- load_nanopore_info(file.path(data_dir, "nanopore", "data", "WIV04.rep1.gap.bed"), nanopore_rep1_fl, nanopore_multiple_rep1, "Rep1")
nanopore_JS_rep1 <- nanopore_JS_rep1[!nanopore_JS_rep1$Name %in% nanopore_multiple_rep1$X4,]
nanopore_JS_rep2 <- load_nanopore_info(file.path(data_dir, "nanopore", "data", "WIV04.rep2.gap.bed"), nanopore_rep2_fl, nanopore_multiple_rep2, "Rep2")
nanopore_JS_rep2 <- nanopore_JS_rep2[!nanopore_JS_rep2$Name %in% nanopore_multiple_rep2$X4,]

nanopore_JS <- rbind(nanopore_JS_rep1, nanopore_JS_rep2)
read_pos_rep1 <- data.frame(
  ReadStart=nanopore_bed_rep1$X2,
  ReadEnd=nanopore_bed_rep1$X3,
  Name=as.character(nanopore_bed_rep1$X4),
  Source="Rep1"
  
)
read_pos_rep2 <- data.frame(
  ReadStart=nanopore_bed_rep2$X2,
  ReadEnd=nanopore_bed_rep2$X3,
  Name=as.character(nanopore_bed_rep2$X4),
  Source="Rep2"
)
read_pos_rep1 <- read_pos_rep1[! read_pos_rep1$Name %in% nanopore_multiple_rep1,]
read_pos_rep2 <- read_pos_rep2[! read_pos_rep2$Name %in% nanopore_multiple_rep2,]
read_pos <- rbind(read_pos_rep1, read_pos_rep2)

## Filter Nanopore JS by position (no care about #Read)
ORF1ab_end <- 21555
leader_orf1ab_df <- nanopore_JS[(nanopore_JS$Start<ORF1ab_end) & (nanopore_JS$End>ORF1ab_end),]
leader_orf1ab_df <- inner_join(leader_orf1ab_df, read_pos)

min_overhang_len <- 20
leader_orf1ab_df <- leader_orf1ab_df[
  ((leader_orf1ab_df$Start-leader_orf1ab_df$ReadStart)>=min_overhang_len) &
  ((leader_orf1ab_df$ReadEnd - leader_orf1ab_df$End)>=min_overhang_len),
  ]
leader_orf1ab_df$GapType <- "LeaderSequence"
leader_orf1ab_df$GapType[leader_orf1ab_df$Start>100] <- "ORF1ab_S2N"

write_tsv(leader_orf1ab_df, file.path(out_dir, "long_JS", "Nanopore.long_JS.sourceData.tsv"))

nanopore_bed_rep1_ORF1ab <- nanopore_bed_rep1[nanopore_bed_rep1$X4 %in% leader_orf1ab_df$Name[leader_orf1ab_df$GapType=="ORF1ab_S2N" & leader_orf1ab_df$Source=="Rep1"],]
nanopore_bed_rep2_ORF1ab <- nanopore_bed_rep2[nanopore_bed_rep2$X4 %in% leader_orf1ab_df$Name[leader_orf1ab_df$GapType=="ORF1ab_S2N" & leader_orf1ab_df$Source=="Rep2"],]
nanopore_bed_merge_ORF1ab <- rbind(nanopore_bed_rep1_ORF1ab, nanopore_bed_rep2_ORF1ab)
write_tsv(nanopore_bed_rep1_ORF1ab, file.path(out_dir, "long_JS", "ORF1ab.rep1.bed12"), col_names = FALSE)
write_tsv(nanopore_bed_rep2_ORF1ab, file.path(out_dir, "long_JS", "ORF1ab.rep2.bed12"), col_names = FALSE)
write_tsv(nanopore_bed_merge_ORF1ab, file.path(out_dir, "long_JS", "ORF1ab.all.bed"), col_names = FALSE)

leader_orf1ab_df$Degraded_5 <- leader_orf1ab_df$ReadStart > fl_cutoff
GapType_info <- leader_orf1ab_df %>% 
  group_by(GapType, Start, End, Degraded_5) %>% 
  summarise(ReadPerEvent=n()) %>%
  group_by(GapType) %>% 
  summarise(
    JSNum=n(), 
    ReadNum=sum(ReadPerEvent), 
    Degrade_5_JSNum=sum(Degraded_5), 
    Degrade_5_ReadNum=sum(ReadPerEvent[Degraded_5])
    )

write_tsv(GapType_info, file.path(out_dir, "long_JS", "Nanopore.long_JS.stat.sourceData.tsv"))
melt_GapType_info <- melt(GapType_info, id.vars = c("GapType"), variable.name = "Stat", value.name = "Number")
melt_GapType_info$Stat <- factor(
  melt_GapType_info$Stat, 
  levels = c("JSNum", "ReadNum", "Degrade_5_JSNum", "Degrade_5_ReadNum"),
  labels = c("#Nanopore JS", "#Nanopore read", "#5' degraded JS", "#5' degraded read")
  )

p <- ggplot(melt_GapType_info, aes(x=GapType, y=Number, fill=GapType)) +
  geom_bar(stat="identity", position="dodge", width = 0.7) +
  geom_text(aes(label=Number), size=1.8, hjust=1.1, color="white") +
  theme_bw() +
  scale_fill_brewer(palette="Set1") +
  facet_wrap(~Stat, scale="free_x", nrow = 1) +
  coord_flip() +
  # scale_y_log10() +
  theme(text = element_text(family="ArialMT", size=6),
        title = element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        panel.grid = element_blank()
  )
ggsave(file.path(out_dir, "long_JS", "Nanopore.long_JS.cnt.pdf"), p, width = 12, height = 2, units = "cm")

melt_leader_orf1ab_df <- melt(
  leader_orf1ab_df[,c("GapType", "Start", "End", "ReadStart", "ReadEnd")], 
  id.vars = c("GapType"), 
  variable.name = "Stat", 
  value.name = "Position"
  )
melt_leader_orf1ab_df$Stat <- factor(
  melt_leader_orf1ab_df$Stat,
  levels = c("Start", "End", "ReadStart", "ReadEnd"),
  labels = c("JS start", "JS end", "Read start", "Read end")
  )

melt_leader_orf1ab_df_start <- melt_leader_orf1ab_df[melt_leader_orf1ab_df$Stat %in% c("JS start", "Read start"), ]
melt_leader_orf1ab_df_end <- melt_leader_orf1ab_df[melt_leader_orf1ab_df$Stat %in% c("JS end", "Read end"), ]

p <- ggplot(melt_leader_orf1ab_df_start, aes(x=Position, fill=GapType)) +
  geom_density(size=0.1) +
  theme_bw() +
  scale_fill_brewer(palette="Set1") +
  facet_wrap(GapType~Stat, scale="free_y", nrow=2) +
  # scale_y_log10() +
  theme(text = element_text(family="ArialMT", size=6),
        title = element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        panel.grid = element_blank()
  )
ggsave(file.path(out_dir, "long_JS", "Nanopore.long_JS.density.start.pdf"), p, width = 6, height = 8, units = "cm")

p <- ggplot(melt_leader_orf1ab_df_end, aes(x=Position, fill=GapType)) +
  geom_density(size=0.1) +
  theme_bw() +
  scale_fill_brewer(palette="Set1") +
  facet_wrap(GapType~Stat, scale="free_y", nrow=2) +
  # scale_y_log10() +
  theme(text = element_text(family="ArialMT", size=6),
        title = element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        panel.grid = element_blank()
  )
ggsave(file.path(out_dir, "long_JS", "Nanopore.long_JS.density.end.pdf"), p, width = 6, height = 8, units = "cm")

# gr <- GRanges("MN996528", IRanges(0, 30000), strand="+")
# orf1ab_bed <- importScore(file.path(out_dir, "long_JS", "ORF1ab.all.bed"), format="BED", ranges=gr)
# orf1ab_bed <- import(file.path(out_dir, "long_JS", "ORF1ab.all.bed"), "bed")
# 
# orf1ab_txdb <- GenomicFeatures::makeTxDbFromGFF(file.path(out_dir, "long_JS", "ORF1ab.all.raw.gff"))
# 
# orf1ab_track <- geneTrack(orf1ab_bed$name, orf1ab_txdb, "gene")
# 
# inche_cm=2.54
# pdf(file.path(out_dir, "long_JS", "ORF1ab.nanopore.pdf"), width=30/inche_cm, height=100/inche_cm)
# viewTracks(trackList(orf1ab_track), gr=gr)
# dev.off()

orf1ab_bed <- import(file.path(out_dir, "long_JS", "ORF1ab.all.bed"), "bed")
orf1ab_bed <- orf1ab_bed[order(orf1ab_bed@ranges@start),]
max_window <- 30000
window_size <- 50

is_hit <- function(start, blocks, width){
  x <- coverage(blocks, shift=(-1*start), width = width)
  res <- (sum(x@values * x@lengths) / width) > 0
  return(res)
}

build_matrix <- function(indx, bed, max_window=30000, window_size=50){
  rec <- bed[indx,]
  start_pos <- start(rec)
  start_li <- seq(0, max_window, window_size)
  hit_var <- sapply(start_li-start_pos, is_hit, blocks=rec$blocks[[1]], width=window_size)
  return(hit_var)
}

hit_mat <- t(sapply(1:length(orf1ab_bed), build_matrix, bed=orf1ab_bed, max_window=max_window, window_size=window_size))
hit_mat <- t(apply(hit_mat, 1, as.character))
ht <- Heatmap(
  hit_mat, 
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_heatmap_legend = FALSE,
  col = c("TRUE"="black", "FALSE"="white"),
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 5),
    title_gp = gpar(fontsize = 5)
  ))

inche_cm=2.54
pdf(file.path(out_dir, "long_JS", "ORF1ab.nanopore.ht.pdf"), width=12/inche_cm, height=20/inche_cm)
print(ht)
dev.off()

inche_cm=2.54
tiff(file.path(out_dir, "long_JS", "ORF1ab.nanopore.ht.tiff"), width=12, height=20, units = "cm", res=600)
print(ht)
dev.off()


