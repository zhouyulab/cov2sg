library(readr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(ComplexHeatmap)

data_dir <- file.path("data", "subgenome")
out_dir <- file.path("analysis", "subgenome")
if (!dir.exists(out_dir)) dir.create(out_dir)

################################################

## Figure 1B

polyA_tail_rep1 <- read_delim("data/polyA_tail/Nanopore_rep1.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
polyA_tail_rep2 <- read_delim("data/polyA_tail/Nanopore_rep2.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
polyA_tail <- rbind(polyA_tail_rep1, polyA_tail_rep2)
polyA_tail$Sample <- "VERO"
polyA_tail$TP <- "48h"

for(tp in c("6h", "12h", "24h")){
  tmp_polyA_tail <- read_delim(sprintf("data/polyA_tail/VERO_%s.tsv", tp), "\t", escape_double = FALSE, trim_ws = TRUE)
  tmp_polyA_tail$Sample <- "VERO"
  tmp_polyA_tail$TP <- tp
  polyA_tail <- rbind(polyA_tail, tmp_polyA_tail)
}

for(tp in c("12h", "24h")){
  tmp_polyA_tail <- read_delim(sprintf("data/polyA_tail/Caco2_%s.tsv", tp), "\t", escape_double = FALSE, trim_ws = TRUE)
  tmp_polyA_tail$Sample <- "Caco2"
  tmp_polyA_tail$TP <- tp
  polyA_tail <- rbind(polyA_tail, tmp_polyA_tail)
}

polyA_tail$Sample <- factor(polyA_tail$Sample, levels = c("VERO", "Caco2"))
polyA_tail$TP <- factor(polyA_tail$TP, levels = c("6h", "12h", "24h", "48h"))
p <- ggplot(polyA_tail, aes(x=Length, color=TP, lty=Sample)) +
  geom_density(fill=NA) +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  labs(x="polyA-tail length", y="#Reads") +
  lims(x=c(0, 200)) +
  theme(text = element_text(family="ArialMT", size=6),
        title = element_text(family="ArialMT", size=6),
        axis.text = element_text(color = "black"),
        panel.grid = element_blank()
  )
ggsave(file.path(out_dir, "seq_pair", "polyA.len.pdf"), p, width = 12, height = 6, units = "cm")


################################################

## Figure 3E


Pair_info_all <- read_delim("analysis/subgenome/seq_pair/Pair.info.all.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
plot_pair_bin <- function(df, bin_size=20){
  df$StartBin <- as.integer(df$Start/bin_size)
  df$EndBin <- as.integer(df$End/bin_size)
  
  bin_info <- df %>% group_by(PairType, StartBin, EndBin) %>% summarise(ReadNum=sum(ReadNum))
  bin_info$Start <- (bin_info$StartBin+0.5) * bin_size
  bin_info$End <- (bin_info$EndBin+0.5) * bin_size
  bin_info <- bin_info[order(bin_info$ReadNum),]
  DL_AL <- bin_info[bin_info$PairType=="DL-AL",]
  DR_AR <- bin_info[bin_info$PairType=="DR-AR",]
  
  p <- ggplot() +
    geom_curve(DR_AR, mapping = aes(x=Start, y=0, xend=End, yend=0, color=log10(ReadNum)), curvature=0.5, size=0.1) +
    geom_curve(DL_AL, mapping = aes(x=Start, y=0, xend=End, yend=0, color=log10(ReadNum)), curvature=-0.5, size=0.1) +
    theme_bw() +
    scale_color_gradient(low="blue", high="red") +
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
consistent_df <- Pair_info_all[Pair_info_all$IsConsistent,]

p <- plot_pair_bin(consistent_df, 20)
ggsave(file.path(out_dir, "seq_pair", "JS.gap_type.curve.consistent.bin_20.pdf"), p, width = 12, height = 6, units = "cm")


################################################

## Figure S6C


all_consistent_gap <- read_delim(file.path(out_dir, "NGS_Nanopore", "all_consistent_gap.tsv"), 
                                 "\t", escape_double = FALSE, trim_ws = TRUE)

NGS_nanpore_consistent_cnt <- left_join(consistent_df, all_consistent_gap)
NGS_nanpore_consistent_cnt$Nanopore_cnt <- NGS_nanpore_consistent_cnt$Nanopore_rep1_Cnt + NGS_nanpore_consistent_cnt$Nanopore_rep2_Cnt
NGS_nanpore_consistent_group_cnt <- NGS_nanpore_consistent_cnt %>% group_by(Group, PairType) %>%
  summarise(ReadNum=sum(ReadNum), Nanopore_cnt=sum(Nanopore_cnt))

cor_df <- NGS_nanpore_consistent_group_cnt %>% group_by(PairType) %>% summarise(cor=cor(ReadNum, Nanopore_cnt, method = "spearman"))
cor_df$Label <- sprintf("R=%.2f", cor_df$cor)
p <- ggplot(NGS_nanpore_consistent_group_cnt, aes(y=log10(ReadNum), x=log10(Nanopore_cnt))) +
  geom_smooth(method = "lm", size=0.3) +
  geom_point(size=0.3) +
  geom_text(data=cor_df, aes(x=1, y=6, label=Label), size=1.8, color="black") +
  scale_x_continuous(breaks = c(1, 3, 5), labels = c("1E1", "1E3", "1E5")) +
  scale_y_continuous(breaks = c(2, 4, 6), labels = c("1E2", "1E4", "1E6")) +
  facet_grid(~PairType) +
  theme_bw() +
  labs(x="#Nanopore reads", y="#NGS reads") +
  theme(text = element_text(family="ArialMT", size=6),
        title = element_text(family="ArialMT", size=6),
        axis.text = element_text(color = "black"),
        legend.key.size = unit(4, "mm"),
        # legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(family="ArialMT", size=6),
        panel.grid = element_blank()
  )
ggsave(file.path(out_dir, "seq_pair", "JS.gap_type.cnt.consistent.pdf"), p, width = 10, height = 4, units = "cm")


################################################

## Figure 3F

junc_cnt_tp <- read_delim("analysis/subgenome/timepoint/JuncNum.SourceData.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
junc_cnt_tp <- junc_cnt_tp[junc_cnt_tp$Sample %in% c("6h", "12h", "24h"),]

NarryKim_df <- read_delim("analysis/subgenome/other_data/Nanopore/NarryKim/cutoff.VeroInf24h.sourceData.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
NarryKim_df <- NarryKim_df[NarryKim_df$FL=="FL", c("Start", "End", "ReadNum", "ReadRatio")]
NarryKim_df$Sample <- "Narry Kim"
NarryKim_df <- NarryKim_df[,c("Sample", "Start", "End", "ReadNum", "ReadRatio")]

caco2_df <- read_delim("analysis/subgenome/caco2/JuncNum.SourceData.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
caco2_df$TP <- sprintf("caco2 %s", caco2_df$TP)
names(caco2_df)[1] <- "Sample"

vero_48h_df <- read_delim("analysis/subgenome/NGS_Nanopore/all_consistent_gap.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
vero_48h_df$Sample <- "48h"
vero_48h_df$ReadNum <- vero_48h_df$Nanopore_rep1_Cnt + vero_48h_df$Nanopore_rep2_Cnt
vero_48h_df$ReadRatio <- vero_48h_df$ReadNum / sum(vero_48h_df$ReadNum)
vero_48h_df <- vero_48h_df[,names(junc_cnt_tp)]

other_js_cnt <- rbind(junc_cnt_tp, vero_48h_df, NarryKim_df, caco2_df)

DL_AL_consistent <- consistent_df[consistent_df$PairType=="DL-AL",c("Start", "End", "ReadNum")]
names(DL_AL_consistent)[3] <- "OriReadNum"
junc_cnt_tp_DL_AL <- left_join(DL_AL_consistent, other_js_cnt)
junc_cnt_tp_DL_AL$ReadRatio <- NULL
junc_cnt_tp_DL_AL$Sample <- factor(junc_cnt_tp_DL_AL$Sample, levels = c("6h", "12h", "24h", "48h", "Narry Kim", "caco2 12h", "caco2 24h"))
junc_cnt_tp_DL_AL_total <- junc_cnt_tp_DL_AL %>% group_by(Start, End) %>% summarise(ReadNum=sum(ReadNum))
junc_cnt_tp_DL_AL_total <- left_join(DL_AL_consistent, junc_cnt_tp_DL_AL_total)
junc_cnt_tp_DL_AL_total <- junc_cnt_tp_DL_AL_total[order(junc_cnt_tp_DL_AL_total$ReadNum, decreasing = T),]
junc_cnt_tp_DL_AL_total$Name <- sprintf("%s-%s", junc_cnt_tp_DL_AL_total$Start, junc_cnt_tp_DL_AL_total$End)

tp_DL_AL_cnt_Mat_df <- dcast(junc_cnt_tp_DL_AL, Start+End+OriReadNum~Sample, value.var = "ReadNum")
tp_DL_AL_cnt_Mat_df <- left_join(DL_AL_consistent, tp_DL_AL_cnt_Mat_df)
tp_DL_AL_cnt_Mat_df <- tp_DL_AL_cnt_Mat_df[order(tp_DL_AL_cnt_Mat_df$OriReadNum, decreasing = T),]
tp_DL_AL_cnt_Mat_df$`NA` <- NULL
tp_DL_AL_hit_mat <- as.matrix(tp_DL_AL_cnt_Mat_df[,4:ncol(tp_DL_AL_cnt_Mat_df)])
tp_DL_AL_hit_mat[is.na(tp_DL_AL_hit_mat)] <- 0
tp_DL_AL_hit_mat <- apply(tp_DL_AL_hit_mat>0, 1, as.character)
colnames(tp_DL_AL_hit_mat) <- sprintf("%s-%s", tp_DL_AL_cnt_Mat_df$Start, tp_DL_AL_cnt_Mat_df$End)
rownames(tp_DL_AL_hit_mat) <- names(tp_DL_AL_cnt_Mat_df[4:ncol(tp_DL_AL_cnt_Mat_df)])

indx_df <- data.frame(Name=junc_cnt_tp_DL_AL_total$Name, RankIndx=1:length(junc_cnt_tp_DL_AL_total$Name))
mat_indx_df <- data.frame(Name=colnames(tp_DL_AL_hit_mat), MatIndx=1:ncol(tp_DL_AL_hit_mat))
indx_df <- left_join(indx_df, mat_indx_df)
tp_DL_AL_hit_mat <- tp_DL_AL_hit_mat[,indx_df$MatIndx]
colnames(tp_DL_AL_hit_mat) == indx_df$Name

tp_libSize_df <- read_delim("analysis/subgenome/timepoint/ReadNum.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
tp_libSize_df <- rbind(
  tp_libSize_df, 
  data.frame(
    Sample="48h", 
    TotalReads=sum(tp_libSize_df$TotalReads[1:2]), 
    FullLengthReads=sum(tp_libSize_df$FullLengthReads[1:2]), 
    FullLengthRatio=sum(tp_libSize_df$FullLengthRatio[1:2])
    ))
tp_libSize_df <- inner_join(data.frame(Sample=rownames(tp_DL_AL_hit_mat)), tp_libSize_df)
narrykim_FL_num <- 229412
caco2_libSize_df <- read_delim("analysis/subgenome/caco2/ReadNum.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)


raw_ha <- rowAnnotation(
  LibSize = anno_barplot(
    c(tp_libSize_df$FullLengthReads, narrykim_FL_num, caco2_libSize_df$FullLengthReads), 
    width = unit(0.8, "cm"), 
    annotation_name_gp=gpar(fontsize = 5),
    axis_param = list(direction = "reverse")
  )
)
col_ha <- columnAnnotation(
  ReadNum=anno_barplot(
    log10(junc_cnt_tp_DL_AL_total$ReadNum), 
    width = unit(0.8, "cm"), 
    annotation_name_gp=gpar(fontsize = 5))
)
ht <- Heatmap(
  tp_DL_AL_hit_mat, 
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  # show_row_names = FALSE,
  show_column_names = FALSE,
  show_heatmap_legend = FALSE,
  col = c("FALSE"="white", "TRUE"="black"),
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 5),
    title_gp  = gpar(fontsize = 5),
    text_gp   = gpar(fontsize = 5)
  ),
  top_annotation = col_ha,
  left_annotation = raw_ha
)

inche_cm <- 2.54
pdf(file.path(out_dir, "seq_pair", "DLAL.consistent.hit.pdf"), width=10/inche_cm, height=5/inche_cm)
print(ht)
dev.off()

################################################


## Figure 1I
JS_cnt_df <- read_delim("analysis/subgenome/timepoint/consistent.SourceData.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
JS_cnt_df$`48h` <- JS_cnt_df$Nanopore_rep1_Cnt + JS_cnt_df$Nanopore_rep2_Cnt
JS_cnt_df <- JS_cnt_df[,c("Start", "End", "6h", "12h", "24h", "48h")]

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

leader_js_cnt_df <- leader_js_df[,c("SubgenomeType", "6h", "12h", "24h", "48h")]
melt_leader_js_cnt_df <- melt(leader_js_cnt_df, id.vars = "SubgenomeType", variable.name = "Time", value.name = "ReadNum")
melt_leader_js_cnt_df$ReadNum <- as.integer(melt_leader_js_cnt_df$ReadNum)
library_size <- melt_leader_js_cnt_df %>% 
  group_by(Time) %>% 
  summarise(TotalRead=sum(ReadNum, na.rm = TRUE))
subgenome_type_df <- melt_leader_js_cnt_df %>% 
  group_by(SubgenomeType, Time) %>% 
  summarise(ReadNum=sum(ReadNum, na.rm = TRUE))
subgenome_type_df$SubgenomeType <- factor(subgenome_type_df$SubgenomeType, levels = c("Other", gene_ref_df$Name[gene_ref_df$Name!="orf1ab"]))
subgenome_type_df$Time <- factor(subgenome_type_df$Time, levels = c("6h", "12h", "24h", "48h"))
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
ggsave(file.path("analysis/subgenome/timepoint", "subgenome.type.cnt.tp.pdf"), p, width = 5, height = 4, units = "cm")

NGS_JS_cnt_df <- read_delim("analysis/subgenome/timepoint/consistent.RNA_seq.SourceData.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
all_consistent_gap <- read_delim("analysis/subgenome/NGS_Nanopore/all_consistent_gap.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
NGS_JS_cnt_df$`48h` <- all_consistent_gap$NGS_rep1_Cnt + all_consistent_gap$NGS_rep2_Cnt
NGS_JS_cnt_df <- NGS_JS_cnt_df[,c("Start", "End", "6h", "12h", "24h", "48h")]
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

NGS_leader_js_cnt_df <- NGS_leader_js_df[,c("SubgenomeType", "6h", "12h", "24h", "48h")]
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
ggsave(file.path("analysis/subgenome/timepoint", "subgenome.type.cnt.tp.NGS.pdf"), p, width = 5, height = 4, units = "cm")

################################################
## Figure S8G-H
Ribo_SARS2_stat <- read_delim("E:/nCoV_VERO/data/Ribo.SARS2.stat.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
Ribo_SARS2_stat$Type <- "Other"
Ribo_SARS2_stat$Type[Ribo_SARS2_stat$`5' site`<100 ] <- "Leader"
Ribo_SARS2_stat$Type[Ribo_SARS2_stat$`5' site`>100 & Ribo_SARS2_stat$`5' site`<20000 & Ribo_SARS2_stat$`3' site` > 20000] <- "ORF1ab"
Ribo_SARS2_stat$Type[Ribo_SARS2_stat$`5' site`>20000 & Ribo_SARS2_stat$`5' site`>20000] <- "S-N"
Ribo_SARS2_info <- Ribo_SARS2_stat %>% group_by(Type) %>% summarise(Num=n())
Ribo_SARS2_info$Type <- factor(Ribo_SARS2_info$Type, levels = c("Leader", "ORF1ab", "S-N"))
Ribo_SARS2_info$y <- -1*sum(Ribo_SARS2_info$Num) + cumsum(Ribo_SARS2_info$Num) - Ribo_SARS2_info$Num / 2
Ribo_SARS2_info$Label <- sprintf("%s\nn=%d", Ribo_SARS2_info$Type, Ribo_SARS2_info$Num)
p <- ggplot(Ribo_SARS2_info, aes(x="", y=-Num, fill=Type)) +
  theme_bw() +
  geom_bar(stat = "identity", width = 1) +
  geom_text(mapping = aes(x=1, y=y, label=Label), size=1.8, color="black") +
  coord_polar("y", start = 0) +
  annotate("text", 1.7, -sum(Ribo_SARS2_stat$Num)/2, label = " ", size=2) +
  scale_fill_brewer(palette = "Set1")+
  theme(
    text = element_text(size = 6),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.spacing = unit(0, units = "mm"),
    plot.margin = margin(t=0.2, r=0, b=0.5, l=2, unit = "cm")
  )
p
ggsave("E:/nCoV_VERO/analysis/subgenome/Ribo/Ribo.type.num.pdf", p, width = 10, height = 5, units = "cm")

Ribo_ORF1 <- Ribo_SARS2_stat[Ribo_SARS2_stat$Type=="ORF1ab",]
Ribo_SN <- Ribo_SARS2_stat[Ribo_SARS2_stat$Type=="S-N",]
p <- ggplot() +
  geom_curve(Ribo_ORF1, mapping = aes(x=`5' site`, y=0, xend=`3' site`, yend=0), color="black", curvature=0.5, size=0.1) +
  geom_curve(Ribo_SN, mapping = aes(x=`5' site`, y=0, xend=`3' site`, yend=0), color="black", curvature=-0.5, size=0.1) +
  theme_bw() +
  scale_x_continuous(breaks = c(0, 10000, 20000, 30000), labels = c("0", "10k", "20k", "30k"), limits = c(0, 29891)) +
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
ggsave("E:/nCoV_VERO/analysis/subgenome/Ribo/Ribo.type.curve.pdf", p, width = 10, height = 5, units = "cm")