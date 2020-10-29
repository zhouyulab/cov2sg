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

ORF1_df <- leader_orf1ab_df[leader_orf1ab_df$GapType=="ORF1ab_S2N",]
SARS2_prot <- read_delim("data/subgenome/genome/SARS2.prot.bed", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
SARS2_prot <- SARS2_prot[order(SARS2_prot$End),]
SARS2_prot <- SARS2_prot[,2:4]
names(SARS2_prot) <- c("Start", "End", "Name")
SARS2_prot <- SARS2_prot[! SARS2_prot$Name %in% c("pp1ab", "pp1a"),]
ORF1_df$ORF1_type <- sapply(ORF1_df$Start, function(x){
  indx <- which.max(which(SARS2_prot$End<x))
  if (length(indx)) {
    res <- SARS2_prot$Name[indx]
  }else{
    res <- "No protein"
  }
  return(res)
})
write_tsv(ORF1_df, file.path(out_dir, "long_JS", "Nanopore.ORF1_JS.sourceData.tsv"))

################################################
## Figure 6E

orf1_junc_num <- ORF1_df %>% group_by(ORF1_type) %>% summarise(ReadNum=n())
SARS2_prot <- SARS2_prot[order(SARS2_prot$End),]
prot_order <- data.frame(ORF1_type=SARS2_prot$Name, Order=1:nrow(SARS2_prot))
orf1_junc_num <- left_join(orf1_junc_num, prot_order)
orf1_junc_num <- orf1_junc_num[order(orf1_junc_num$Order, decreasing = T),]
orf1_junc_num$CumSumReadNum <- cumsum(orf1_junc_num$ReadNum)
orf1_junc_num$CumSumReadNum[orf1_junc_num$ORF1_type=="No protein"] <- NA
orf1_junc_num$Order[orf1_junc_num$ORF1_type=="No protein"] <- max(orf1_junc_num$Order, na.rm = T) + 1
orf1_junc_num <- orf1_junc_num[order(orf1_junc_num$Order, decreasing = T),]
orf1_junc_num$Name <- factor(1:nrow(orf1_junc_num), levels = 1:nrow(orf1_junc_num), labels = (orf1_junc_num$ORF1_type))
melt_orf1_junc <- melt(orf1_junc_num[,c("Name", "ReadNum", "CumSumReadNum")], id.vars = "Name", variable.name = "Type", value.name = "ReadNum")
melt_orf1_junc$Type <- factor(melt_orf1_junc$Type, levels = c("ReadNum", "CumSumReadNum"), labels = c("#Nanopore read", "Cumulative #read"))
p <- ggplot(melt_orf1_junc, aes(x=Name, y=ReadNum, fill=Name)) +
  geom_bar(stat="identity", position="dodge", width = 0.7) +
  geom_text(aes(label=ReadNum), size=1.8, hjust=1.1) +
  theme_bw() +
  # scale_fill_brewer(palette="Set1") +
  facet_wrap(~Type, scale="free_x") +
  coord_flip() +
  scale_y_log10() +
  theme(text = element_text(family="ArialMT", size=6),
        title = element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        panel.grid = element_blank()
  )
ggsave(file.path(out_dir, "long_JS", "Nanopore.ORF1_JS.num.pdf"), p, width = 5, height = 4.7, units = "cm")


################################################
## Figure 6D

start_cnt <- ORF1_df %>% group_by(ORF1_type, Start) %>% 
  summarise(ReadNum=n()) %>% 
  group_by(ORF1_type) %>% 
  arrange(-ReadNum, .by_group = TRUE) %>%
  mutate(Major=ReadNum==max(ReadNum))
start_cnt <- start_cnt[start_cnt$Start>100,]
y_indx <- data.frame(ORF1_type=c(as.character(orf1_junc_num$Name[orf1_junc_num$Name!="No protein"]), "No protein"), y_indx=1:nrow(orf1_junc_num))

start_cnt <- left_join(start_cnt, y_indx)
start_cnt$y_end <- start_cnt$y_indx + 0.3
start_cnt$y_end[start_cnt$Major] <- start_cnt$y_indx[start_cnt$Major] + 0.5

y_indx$Name <- y_indx$ORF1_type
tmp_SARS2_prot <- SARS2_prot
tmp_SARS2_prot$ORF1_type <- tmp_SARS2_prot$Name
tmp_SARS2_prot <- inner_join(tmp_SARS2_prot, y_indx)
tmp_SARS2_prot$Mid <- (tmp_SARS2_prot$Start + tmp_SARS2_prot$End) / 2
max_site_df <- start_cnt %>% group_by(ORF1_type, y_indx) %>% summarise(max_site=max(Start))

p <- ggplot() +
  geom_linerange(data=max_site_df, mapping = aes(xmin=0, xmax=max_site, y=y_indx), color="black", size=0.3) +
  geom_linerange(data=start_cnt[start_cnt$ReadNum<=10,], mapping = aes(x=Start, ymin=y_end, ymax=y_indx), color="grey70", size=0.3) +
  geom_linerange(data=start_cnt[start_cnt$ReadNum>10,], mapping = aes(x=Start, ymin=y_end, ymax=y_indx), color="red", size=0.3) +
  geom_rect(data=tmp_SARS2_prot, aes(xmin=Start, xmax=End, ymin=y_indx-0.1, ymax=y_indx+0.1), fill="black") +
  geom_text(data=tmp_SARS2_prot, aes(x=Mid, y=y_indx+0.2, label=Name), color="black", size=1.2) +
  theme_bw() +
  theme(text = element_text(family="ArialMT", size=5),
        title = element_text(family="ArialMT", size=5),
        axis.text = element_blank(),
        panel.border = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()
  )
p
ggsave(file.path(out_dir, "long_JS", "ORF1.start.pos.pdf"), p, width = 10, height = 8, units = "cm")


nanopore_bed_rep1_ORF1ab <- nanopore_bed_rep1[nanopore_bed_rep1$X4 %in% leader_orf1ab_df$Name[leader_orf1ab_df$GapType=="ORF1ab_S2N" & leader_orf1ab_df$Source=="Rep1"],]
nanopore_bed_rep2_ORF1ab <- nanopore_bed_rep2[nanopore_bed_rep2$X4 %in% leader_orf1ab_df$Name[leader_orf1ab_df$GapType=="ORF1ab_S2N" & leader_orf1ab_df$Source=="Rep2"],]

write_tsv(nanopore_bed_rep1_ORF1ab, file.path(out_dir, "long_JS", "ORF1ab.rep1.bed12"), col_names = FALSE)
write_tsv(nanopore_bed_rep2_ORF1ab, file.path(out_dir, "long_JS", "ORF1ab.rep2.bed12"), col_names = FALSE)

################################################
## Figure S8F

leader_orf1ab_df$Degraded_5 <- leader_orf1ab_df$ReadStart > fl_cutoff
GapType_info <- leader_orf1ab_df %>% 
  group_by(GapType, Start, End, Degraded_5) %>% 
  summarise(ReadPerEvent=n()) %>%
  group_by(GapType) %>% 
  summarise(
    JSNum=n(), 
    ReadNum=sum(ReadPerEvent)
    )

write_tsv(GapType_info, file.path(out_dir, "long_JS", "Nanopore.long_JS.stat.sourceData.tsv"))
melt_GapType_info <- melt(GapType_info, id.vars = c("GapType"), variable.name = "Stat", value.name = "Number")
melt_GapType_info$Stat <- factor(
  melt_GapType_info$Stat, 
  levels = c("JSNum", "ReadNum"),
  labels = c("#Nanopore JS", "#Nanopore read")
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

################################################
## Figure S8E

melt_leader_orf1ab_df <- melt(
  leader_orf1ab_df[,c("GapType", "Start", "End")], 
  id.vars = c("GapType"), 
  variable.name = "Stat", 
  value.name = "Position"
  )
melt_leader_orf1ab_df$Stat <- factor(
  melt_leader_orf1ab_df$Stat,
  levels = c("Start", "End"),
  labels = c("JS start", "JS end")
  )

melt_leader_orf1ab_df_start <- melt_leader_orf1ab_df[melt_leader_orf1ab_df$Stat %in% c("JS start"), ]
melt_leader_orf1ab_df_end <- melt_leader_orf1ab_df[melt_leader_orf1ab_df$Stat %in% c("JS end"), ]

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


################################################
## Figure 5C


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


