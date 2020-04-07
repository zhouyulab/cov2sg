library(readr)
library(dplyr)
library(ggplot2)
library(BiocParallel)
library(ComplexHeatmap)
library(RColorBrewer)
library(reshape2)
library(igraph)
library(circlize)


data_dir <- file.path("data", "subgenome")
out_dir <- file.path("analysis", "subgenome")
if (!dir.exists(out_dir)) dir.create(out_dir)

#' Read RNAcofold predicted MFE structure between fusion segments
load_pair <- function(f_loc, f_pair, NGS_JS=NULL){
  seq_info <- read_delim(f_loc, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  seq_info <- seq_info[, c(2, 3, 4, 8, 9, 10)]
  names(seq_info) <- c("Start", "End", "Name", "ReadNum", "LeftSeq", "RightSeq")
  seq_pair <- read_table2(f_pair, col_names = FALSE)
  seq_pair <- seq_pair[, c(2, 4, 5)]
  seq_pair$X4 <- gsub("\"", "", seq_pair$X4)
  names(seq_pair) <- c("Name", "Pair", "Energy")
  seq_pair$LeftPair <- sapply(strsplit(seq_pair$Pair, "&"), function(x) {return(x[1])})
  seq_pair$RightPair <- sapply(strsplit(seq_pair$Pair, "&"), function(x) {return(x[2])})
  seq_pair <- seq_pair[, c("Name", "LeftPair", "RightPair", "Energy")]
  pair_df <- inner_join(seq_info, seq_pair)
  # check if there is internal pairing in either segment
  pair_df$Paired <- (!grepl("\\)", as.character(pair_df$LeftPair))) & (!grepl("\\(", as.character(pair_df$RightPair)))
  if (!is.null(NGS_JS)){
    pair_df <- inner_join(pair_df, NGS_JS)
  }
  pair_df$Name <- NULL
  return(pair_df)
}

#' Merge neighbouring junctions within window into group
creat_group <- function(gap_df, resolution=5){
  gap_df <- gap_df[, c("Start", "End", "ReadNum")]
  gap_df <- gap_df[order(gap_df$Start, gap_df$End), ]
  is_one_group <- function(junc_df1, junc_df2, resolution=5){
    diff_junc <- junc_df1 - junc_df2
    return(all(abs(diff_junc) <= resolution))
  }
  indx1 <- c()
  indx2 <- c()
  for (i in 1:nrow(gap_df)) {
    junc_df1 <- gap_df[i,]
    for (j in i:nrow(gap_df)) {
      junc_df2 <- gap_df[j, ]
      if ((junc_df2$Start-junc_df1$Start) > resolution) break
      can_merge <- is_one_group(junc_df1[, c("Start", "End")], junc_df2[, c("Start", "End")], resolution)
      if(can_merge){
        indx1 <- c(indx1, i)
        indx2 <- c(indx2, j)
      }
    }
  }
  edge_df <- data.frame(from=indx1, to=indx2)
  edge_df <- edge_df[order(edge_df$from, edge_df$to), ]
  subgenome_graph <- graph_from_edgelist(as.matrix(edge_df), directed = FALSE)
  graph_cluster <- components(subgenome_graph, "strong")
  
  merge_gap <- function(indx_li, gap_df){
    selected_df <- gap_df[indx_li, ]
    selected_df <- selected_df[order(selected_df$ReadNum, decreasing = TRUE), ]
    res <- data.frame(
      Start = selected_df$Start,
      End = selected_df$End,
      Group = sprintf("%d-%d", selected_df$Start[1], selected_df$End[1]), # Max ReadNum as core
      LeaderStart = selected_df$Start[1],
      LeaderEnd = selected_df$End[1],
      GroupSize = nrow(selected_df),
      GroupReadNum = sum(selected_df$ReadNum),
      ReadNum = selected_df$ReadNum
    )
    return(res)
  }
  merged_group <- do.call(rbind, lapply(groups(graph_cluster), merge_gap, gap_df=gap_df))
  merged_group <- merged_group[order(merged_group$GroupReadNum, decreasing = TRUE), ]
  return(merged_group)
}

#' Read junction sites
load_NGS_js <- function(f){
  NGS_JS <- read_delim(f, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  NGS_JS <- NGS_JS[, c(4, 5, 8)]
  names(NGS_JS) <- c("Start", "End", "Strand")
  NGS_JS <- NGS_JS[NGS_JS$Strand == "+", ]
  NGS_JS$Strand <- NULL
  return(NGS_JS)
}

##     Donor site    Acceptor site
##      XXXDL|DR........AL|ARYYY

## Forward (DL-AL pair)
## From        XXXDL|
## To             AL|ARYYY
## Subgenome:  XXXDL|DRYYY


## Reverse (DR-AR pair)
## From             |ARYYY
## To          XXXDL|DR
## Subgenome:  XXXDL|ARYYY

js_data_dir <- file.path(data_dir, "js")
NGS_JS_rep1 <- load_NGS_js(file.path(js_data_dir, "CoV_rep1_jd.csv"))
NGS_JS_rep2 <- load_NGS_js(file.path(js_data_dir, "CoV_rep2_jd.csv"))
both_NGS_JS <- inner_join(NGS_JS_rep1, NGS_JS_rep2)

l_df <- load_pair(file.path(data_dir, "seq_pair/CoV_merged_lft.csv"), 
                  file.path(data_dir, "seq_pair/CoV_merged_lft_hybcons.csv"), both_NGS_JS)
r_df <- load_pair(file.path(data_dir, "seq_pair/CoV_merged_rgt.csv"),
                  file.path(data_dir, "seq_pair/CoV_merged_rgt_hybcons.csv"), both_NGS_JS)

df <- inner_join(l_df[, c("Start", "End")], r_df[, c("Start", "End")])
l_df <- left_join(df, l_df)
r_df <- left_join(df, r_df)
names(l_df)[4:9] <- c("DL_seq", "AL_seq", "DL_AL_pair_DL", "DL_AL_pair_AL", "DL_AL_pair_Energy", "has_DL_AL_pair")
names(r_df)[4:9] <- c("DR_seq", "AR_seq", "DR_AR_pair_DR", "DR_AR_pair_AR", "DR_AR_pair_Energy", "has_DR_AR_pair")


# Subgroup pairings by different types
LEADER_BOUNDARY <- 100
ORF1ab_BOUNDARY <- 20000
group_df <- creat_group(l_df, resolution=5)
group_df$GapType <- "Other"
group_df$GapType[group_df$Start <= LEADER_BOUNDARY] <- "LeaderSequence"
group_df$GapType[(group_df$Start > LEADER_BOUNDARY) & (group_df$Start <= ORF1ab_BOUNDARY) & (group_df$End > ORF1ab_BOUNDARY)] <- "ORF1ab_S2N"
group_df$GapType[(group_df$Start > ORF1ab_BOUNDARY) & (group_df$End > ORF1ab_BOUNDARY)] <- "S2N_S2N"
group_df$GapType <- factor(group_df$GapType, levels = c("LeaderSequence", "ORF1ab_S2N", "S2N_S2N", "Other"))

# Read gene annotation
gene_ref <- read_delim(file.path(data_dir, "genome/WIV04.bed"),
                       "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
gene_ref_df <- data.frame(Name=gene_ref$X4, Start=gene_ref$X2, End=gene_ref$X3)
gene_ref_df$Mid <- (gene_ref_df$Start + gene_ref_df$End) / 2
gene_ref_df <- gene_ref_df[order(gene_ref_df$Start), ]
gene_ref_df$Name <- as.character(gene_ref_df$Name)

group_df$SubgenomeType <- sapply(group_df$End, function (end_site) {
  # locate the first downstream complete gene of the first junction
  idx <- which(end_site < gene_ref_df$Start)
  if (length(idx) == 0) {
    sgtype <- "Other"
  } else {
    sgtype <- gene_ref_df$Name[min(idx)] 
  }
  return(sgtype)}
)

oSP <- file.path(out_dir, "seq_pair")
if (!dir.exists(oSP)) dir.create(oSP)

geneS_start <- gene_ref_df[gene_ref_df$Name == "S", "Start"]
group_df$SubgenomeType[group_df$Start > geneS_start] <- "Other" # Subgroup of subgenome downstream of S gene
group_df$SubgenomeType <- factor(group_df$SubgenomeType, levels = c(rev(gene_ref_df$Name), "Other"))

group_df <- group_df[order(group_df$GroupReadNum, group_df$Group, group_df$ReadNum, decreasing = TRUE), ]
merge_df <- left_join(group_df, l_df)
merge_df <- left_join(merge_df, r_df)


## AMGAAC motif
TRScore <- read_delim(file.path(data_dir, "jsmotif/TRScore.bed"), "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
names(TRScore) <- c("Chrom", "Start", "End", "Name", "Score", 'Strand')
TRScore <- TRScore[,c("Start", "End", "Name", "Strand")]
TRScore <- TRScore[order(TRScore$Start),]

start_num <- merge_df %>% group_by(Start) %>% summarise(Number=sum(ReadNum))
names(start_num)[1] <- "Position"
end_num <- merge_df %>% group_by(End) %>% summarise(Number=sum(ReadNum))
names(end_num)[1] <- "Position"

build_neighbor_read_num <- function(line, num_li, expand_size=10){
  motif_start <- line[1]
  motif_end <- line[2]
  min_indx <- motif_start - expand_size
  max_indx <- motif_end + expand_size - 1
  res <- c()
  for(i in min_indx:max_indx){
    if(i %in% num_li$Position){
      res <- c(res, num_li$Number[num_li$Position==i])
    }else{
      res <- c(res, 0)
    }
  }
  return(res)
}

start_mat <- t(apply(TRScore[,c("Start", "End")], 1, build_neighbor_read_num, num_li=start_num))
end_mat <- t(apply(TRScore[,c("Start", "End")], 1, build_neighbor_read_num, num_li=end_num))
colnames(start_mat) <- c(seq(-10, -1), "A", "M", "G", "A", "A", "C", seq(1:10))
colnames(end_mat) <- c(seq(-10, -1), "A", "M", "G", "A", "A", "C", seq(1:10))
rownames(start_mat) <- paste(TRScore$Start, TRScore$End, sep = "-")
rownames(end_mat) <- paste(TRScore$Start, TRScore$End, sep = "-")

ht1 <- Heatmap(
  log10(start_mat+1), 
  col=colorRamp2(c(0, 8), c("white", "red")),
  column_title = "#Donor site",
  column_title_gp  = gpar(fontsize = 7),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 5),
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 5),
    title_gp = gpar(fontsize = 5)
  )
)
ht2 <- Heatmap(
  log10(end_mat+1), 
  col=colorRamp2(c(0, 8), c("white", "red")),
  column_title = "#Accetpor site",
  column_title_gp  = gpar(fontsize = 7),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 5),
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 5),
    title_gp = gpar(fontsize = 5)
  )
)
motif_col <- c("#1B9E77", "#D95F02")
names(motif_col) <- c("ACGAAC", "AAGAAC")
strand_col <- c("#7570B3", "#E7298A")
names(strand_col) <- c("+", "-")
ha <- rowAnnotation(
  Motif=TRScore$Name,
  Strand=TRScore$Strand,
  annotation_name_gp=gpar(fontsize = 5),
  width=unit(1, "cm"),
  border=FALSE,
  col=list(
    Motif=motif_col,
    Strand=strand_col
  )
)

inche_cm <- 2.54
pdf(file.path(oSP, "Motif.ht.pdf"), width=14/inche_cm, height=8/inche_cm)
print(ha + ht1 + ht2)
dev.off()

fwd_motif <- TRScore[TRScore$Strand=="+",]
rev_motif <- TRScore[TRScore$Strand=="-",]

has_motif <- function(line, motif_df, expand=10){
  donor_pos <- line[1]
  acceptor_pos <- line[2]
  donor_match <- (motif_df$Start >= (donor_pos-expand)) & (motif_df$End <= (donor_pos+expand))
  acceptor_match <- (motif_df$Start >= (acceptor_pos-expand)) & (motif_df$End <= (acceptor_pos+expand))
  pair_match <- any(donor_match) & any(acceptor_match)
  return(any(pair_match))
}
merge_df$HasMotif <- apply(merge_df[,c("Start", "End")], 1, has_motif, motif_df=fwd_motif) | apply(merge_df[,c("Start", "End")], 1, has_motif, motif_df=rev_motif)
write_tsv(merge_df, file.path(oSP, "seq_pair.SourceData.tsv"))

# Read jumping pairs from random gaps
rand_l_df <- load_pair(file.path(data_dir, "seq_pair/rand_lft.csv"),
                       file.path(data_dir, "seq_pair/rand_lft_hybcons.csv"))
rand_r_df <- load_pair(file.path(data_dir, "seq_pair/rand_rgt.csv"),
                       file.path(data_dir, "seq_pair/rand_rgt_hybcons.csv"))
names(rand_l_df)[4:9] <- c("DL_seq", "AL_seq", "DL_AL_pair_DL", "DL_AL_pair_AL", "DL_AL_pair_Energy", "has_DL_AL_pair")
names(rand_r_df)[4:9] <- c("DR_seq", "AR_seq", "DR_AR_pair_DR", "DR_AR_pair_AR", "DR_AR_pair_Energy", "has_DR_AR_pair")
merge_rand_df <- inner_join(rand_l_df, rand_r_df)
write_tsv(merge_rand_df, file.path(out_dir, "seq_pair", "seq_pair.random.SourceData.tsv"))

# Analyze MFE
merge_df$Name <- sprintf("%d-%d", merge_df$Start, merge_df$End)
energy_cutoff <- 0
s_merge_df <- merge_df[order(merge_df$ReadNum, decreasing = T), ]
s_merge_df <- s_merge_df[s_merge_df$has_DL_AL_pair & s_merge_df$has_DR_AR_pair, ]
s_merge_df$R_strong <- s_merge_df$DL_AL_pair_Energy > s_merge_df$DR_AR_pair_Energy
random_energy_df <- merge_rand_df[merge_rand_df$has_DL_AL_pair & merge_rand_df$has_DR_AR_pair, ]
random_energy_df$R_strong <- random_energy_df$DL_AL_pair_Energy > random_energy_df$DR_AR_pair_Energy

energy_high_df <- melt(s_merge_df[1:200, c("R_strong", "DL_AL_pair_Energy", "DR_AR_pair_Energy")],
                       id="R_strong", variable.name = "Type", value.name = "Energy")
energy_low_df <- melt(s_merge_df[(nrow(s_merge_df)-200):nrow(s_merge_df), c("R_strong", "DL_AL_pair_Energy", "DR_AR_pair_Energy")],
                      id="R_strong", variable.name = "Type", value.name = "Energy")
random_energy_df <- melt(random_energy_df[, c("R_strong", "DL_AL_pair_Energy", "DR_AR_pair_Energy")],
                         id="R_strong", variable.name = "Type", value.name = "Energy")

energy_high_df$Source <- "Junction (High level)"
energy_low_df$Source <- "Junction (Low level)"
random_energy_df$Source <- "Random"

# MFE distribution of Strong/Weak/Random jumping pairs 
merge_energy_df <- rbind(energy_high_df, energy_low_df, random_energy_df)
merge_energy_df$Type <- factor(merge_energy_df$Type, 
                               levels = c("DL_AL_pair_Energy", "DR_AR_pair_Energy"), 
                               labels = c("DL-AL", "DR-AR"))
merge_energy_df$R_strong <- factor(merge_energy_df$R_strong, 
                                   levels = c(TRUE, FALSE), 
                                   labels = c("DR-AR Strong", "DL-AL Strong"))
p <- ggplot(merge_energy_df, aes(x=Energy, color=Source)) +
  geom_density(fill=NA) +
  scale_color_manual(values = c("Junction (High level)"="red", "Junction (Low level)"="blue", "Random"="grey70")) +
  theme_bw() +
  facet_grid(R_strong~Type) +
  theme(text = element_text(family="ArialMT", size=6),
        title = element_text(family="ArialMT", size=6),
        axis.text = element_text(color = "black"),
        axis.ticks.y = element_blank(),
        legend.key.size = unit(4, "mm"),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(family="ArialMT", size=6),
        panel.grid = element_blank()
  )
ggsave(file.path(oSP, "Energy.cutoff.pdf"), p, width = 11, height = 8, units = "cm")


# Convert pairing structure to boolean matrix (T: paired, F: unpaired) 
pair2mat <- function(pair_str_li, name_li, has_paired, energy, energy_cutoff=0){
  mat <- t(sapply(strsplit(pair_str_li, ""), function(x) { return(x) }))
  mat <- mat != "."
  mat[has_paired=="FALSE", ] <- FALSE
  mat[energy>energy_cutoff, ] <- FALSE
  rownames(mat) <- name_li
  return(mat)
}

#' Show pairing pattern
plot_ht <- function(selected_df, energy_cutoff, split_group=NULL, order_indx=NULL, diff_energy_cutoff=2){
  DL_AL_pair_DL_mat <- pair2mat(selected_df$DL_AL_pair_DL, selected_df$Name, selected_df$has_DL_AL_pair, selected_df$DL_AL_pair_Energy, energy_cutoff)
  DL_AL_pair_AL_mat <- pair2mat(selected_df$DL_AL_pair_AL, selected_df$Name, selected_df$has_DL_AL_pair, selected_df$DL_AL_pair_Energy, energy_cutoff)
  DR_AR_pair_DR_mat <- pair2mat(selected_df$DR_AR_pair_DR, selected_df$Name, selected_df$has_DR_AR_pair, selected_df$DR_AR_pair_Energy, energy_cutoff)
  DR_AR_pair_AR_mat <- pair2mat(selected_df$DR_AR_pair_AR, selected_df$Name, selected_df$has_DR_AR_pair, selected_df$DR_AR_pair_Energy, energy_cutoff)
  
  fwd_mat <- cbind(DL_AL_pair_DL_mat, DL_AL_pair_AL_mat)
  rev_mat <- cbind(DR_AR_pair_DR_mat, DR_AR_pair_AR_mat)
  merged_mat <- cbind(DL_AL_pair_DL_mat, DL_AL_pair_AL_mat, DR_AR_pair_DR_mat, DR_AR_pair_AR_mat)
  
  diff_energy <- selected_df$DL_AL_pair_Energy - selected_df$DR_AR_pair_Energy
  blank_condition1 <- apply(merged_mat, 1, sum) == 0  # Unpaired in neither mode
  blank_condition2 <- (abs(diff_energy) <= diff_energy_cutoff) # Not enough MFE difference to classify between two modes
  blank_condition3 <- (apply(DL_AL_pair_DL_mat, 1, sum) == 0 | apply(DL_AL_pair_AL_mat, 1, sum) == 0) & (diff_energy < 0) # Left more stable, but unpaired
  blank_condition4 <- (apply(DR_AR_pair_DR_mat, 1, sum) == 0 | apply(DR_AR_pair_AR_mat, 1, sum) == 0) & (diff_energy > 0) # Right more stable, but unpaired
  blank_indx <- which(blank_condition1 | blank_condition2 | blank_condition3 | blank_condition4) # all are cases without mode assignment
  blank_df <- data.frame(OriIndx = blank_indx)
  
  # clustering if ordering not given 
  if(is.null(order_indx)) {
    ## "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
    ##  "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)
    l_indx <- which(diff_energy < 0)
    r_indx <- which(diff_energy >= 0)
    l_indx <- setdiff(l_indx, blank_indx)
    r_indx <- setdiff(r_indx, blank_indx)
    stopifnot(all(sort(c(l_indx, r_indx, blank_indx)) == 1:nrow(selected_df)))
    
    # Order by read# for undefined mode
    blank_indx_df <- data.frame(OriIndx=blank_indx, ReadNum=selected_df$ReadNum[blank_indx], AllUnpair=blank_condition1[blank_indx])
    blank_indx_df <- blank_indx_df[order(blank_indx_df$AllUnpair, -blank_indx_df$ReadNum), ]
    
    # Define which side to be used in clustering
    l_indx_df <- data.frame(OriIndx=l_indx, SubIndx=1:length(l_indx))
    r_indx_df <- data.frame(OriIndx=r_indx, SubIndx=1:length(r_indx))
    
    l_dist_mat <- dist(fwd_mat[l_indx, ], method="manhattan")
    r_dist_mat <- dist(rev_mat[r_indx, ], method="manhattan")
    l_clu <- hclust(l_dist_mat, method="ward.D2")
    r_clu <- hclust(r_dist_mat, method="ward.D2")
    l_order_df <- data.frame(SubIndx = l_clu$order)
    r_order_df <- data.frame(SubIndx = r_clu$order)
    l_order_df <- left_join(l_order_df, l_indx_df)
    r_order_df <- left_join(r_order_df, r_indx_df)
    # Orderbing by Left, Right, Unknown modes
    order_indx <- c(r_order_df$OriIndx, l_order_df$OriIndx, blank_indx_df$OriIndx)
  }
  merged_row_clu_mat <- merged_mat[order_indx, ]
  merged_row_clu_mat <- t(apply(merged_row_clu_mat, 1, as.character))
  merged_df_row_clu <- selected_df[order_indx, ]
  order_df <- data.frame(NewIndx=1:nrow(selected_df), OriIndx=order_indx)
  blank_df <- inner_join(blank_df, order_df)
  
  #' Find >= 6 continuous paired segments allowing one mismatch
  find_continue_pair <- function(line, kmer=6, mismatch=1){
    line[line=="FALSE"] <- "G"
    line[line=="TRUE"] <- "P"
    cp_indx <- c()
    for (i in 1:length(line)) {
      eidx <- i+kmer-1
      match_indx <- which(line[i:eidx] == "P")
      if (length(match_indx) >= (kmer-mismatch)) {
        cp_indx <- c(cp_indx, i-1+match_indx)
      }
    }
    cp_indx <- unique(cp_indx)
    line[cp_indx] <- "C" # Contunuous paired state
    return(line)
  }
  
  build_ht_mat <- function(mat, order_indx, no_pair){
    new_mat <- t(apply(mat, 1, as.character))
    new_mat <- t(apply(new_mat, 1, find_continue_pair))
    new_mat[no_pair, ] <- "N" # Left and Right both unpaired
    new_mat <- new_mat[order_indx, ]
    return(new_mat)
  }
  no_pair <- rowSums(merged_mat) == 0
  DL_AL_pair_DL_row_clu_mat <- build_ht_mat(DL_AL_pair_DL_mat, order_indx, no_pair)
  DL_AL_pair_AL_row_clu_mat <- build_ht_mat(DL_AL_pair_AL_mat, order_indx, no_pair)
  DR_AR_pair_DR_row_clu_mat <- build_ht_mat(DR_AR_pair_DR_mat, order_indx, no_pair)
  DR_AR_pair_AR_row_clu_mat <- build_ht_mat(DR_AR_pair_AR_mat, order_indx, no_pair)
  
  create_ht_body <- function(mat, title, split_group=NULL){
    raw_ha <- rowAnnotation(
      PairedNumber = anno_barplot(
        apply(mat, 1, function(x) { return(sum(x %in% c("C", "P"))) }), 
        width = unit(0.8, "cm"), 
        annotation_name_gp=gpar(fontsize = 5),
        axis_param = list(direction = "reverse")
        )
    )
    col_ha <- columnAnnotation(
      GapNumber=anno_barplot(
        apply(mat, 2, function(x) { return(sum(x %in% c("C", "P"))) }), 
        width = unit(0.8, "cm"), 
        annotation_name_gp=gpar(fontsize = 5))
    )
    ht <- Heatmap(
      mat, 
      column_title = title,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = FALSE,
      show_column_names = FALSE,
      show_heatmap_legend = FALSE,
      split = split_group,
      col = c("G"="white", "P"="#FEB24C", "C"="#F03B20", "N"="grey70"),
      heatmap_legend_param = list(
        labels_gp = gpar(fontsize = 5),
        title_gp = gpar(fontsize = 5)
      ),
      top_annotation = col_ha,
      left_annotation = raw_ha
    )
    return(ht)
  }
  if (!is.null(split_group)) {
    split_group <- merged_df_row_clu[, split_group]
  }
  ht1 <- create_ht_body(DL_AL_pair_DL_row_clu_mat, "DL-AL pair (DL)", split_group)
  ht2 <- create_ht_body(DL_AL_pair_AL_row_clu_mat, "DL-AL pair (AL)", split_group)
  ht3 <- create_ht_body(DR_AR_pair_DR_row_clu_mat, "DR-AR pair (DR)", split_group)
  ht4 <- create_ht_body(DR_AR_pair_AR_row_clu_mat, "DR-AR pair (AR)", split_group)
  
  diff_energy <- merged_df_row_clu$DL_AL_pair_Energy - merged_df_row_clu$DR_AR_pair_Energy
  diff_energy_tag <- "NC" # Unclassifed
  diff_energy_tag[diff_energy > diff_energy_cutoff] <- "DR-AR"
  diff_energy_tag[diff_energy < (-1*diff_energy_cutoff)] <- "DL-AL"
  diff_energy_tag[blank_df$NewIndx] <- "NC"
  
  GapType_col <- brewer.pal(length(levels(merged_df_row_clu$GapType)), "Set1")
  names(GapType_col) <- levels(merged_df_row_clu$GapType)
  SubgenomeType_col <- c(brewer.pal(length(levels(merged_df_row_clu$SubgenomeType))-1, "RdBu"), "grey70")
  names(SubgenomeType_col) <- levels(merged_df_row_clu$SubgenomeType)
  DiffEnergy_col <- c("DR-AR"="red", "DL-AL"="blue", "NC"="grey70")
  black_bool_col <- c("black", "white")
  names(black_bool_col) <- c("TRUE", "FALSE")
  ha <- rowAnnotation(
    DiffEnergy = diff_energy_tag,
    GapType = merged_df_row_clu$GapType, 
    SubgenomeType = merged_df_row_clu$SubgenomeType,
    Motif = as.character(merged_df_row_clu$HasMotif),
    ForwardEndPair = as.character((DL_AL_pair_DL_row_clu_mat[, ncol(DL_AL_pair_DL_row_clu_mat)] %in% c("C", "P")) &
                                  (DL_AL_pair_AL_row_clu_mat[, 1] %in% c("C", "P"))),
    ReverseEndPair = as.character((DR_AR_pair_DR_row_clu_mat[, 1] %in% c("C", "P")) & 
                                  (DR_AR_pair_AR_row_clu_mat[, ncol(DR_AR_pair_AR_row_clu_mat)] %in% c("C", "P"))),
    LogReadNum = anno_barplot(log10(merged_df_row_clu$ReadNum), width=unit(1.5, "cm"), gp=gpar(fontsize = 5)),
    annotation_name_gp=gpar(fontsize = 5),
    col=list(
      DiffEnergy=DiffEnergy_col,
      GapType=GapType_col,
      SubgenomeType=SubgenomeType_col,
      Motif=black_bool_col,
      ForwardEndPair=black_bool_col,
      ReverseEndPair=black_bool_col
    )
  )
  p <- ht3 + ht4 + ht1 + ht2 + ha
  return(p)
}

# All junctions seen in both NGS Rep1 and Rep2, not care about read number
inche_cm <- 2.54
p_merge <- plot_ht(merge_df, energy_cutoff)
pdf(file.path(oSP, "seq_pair.all.ht.pdf"), width=30/inche_cm, height=(5+2*log2(nrow(merge_df)))/inche_cm)
print(p_merge)
dev.off()
cat("seq_pair.all#:", nrow(merge_df))

p_merge_by_gap <- plot_ht(merge_df, energy_cutoff, split_group="GapType")
pdf(file.path(oSP, "seq_pair.all.by_gap.ht.pdf"), width=30/inche_cm, height=(5+2*log2(nrow(merge_df)))/inche_cm)
print(p_merge_by_gap)
dev.off()

merge_rand_df$GapType <- "Other"
merge_rand_df$SubgenomeType <- "Other"
merge_rand_df$GapType <- factor(merge_rand_df$GapType, levels = levels(merge_df$GapType))
merge_rand_df$SubgenomeType <- factor(merge_rand_df$SubgenomeType, levels = levels(merge_df$SubgenomeType))
merge_rand_df$Name <- sprintf("%s-%s", merge_rand_df$Start, merge_rand_df$End)
merge_rand_df$ReadNum <- 2
merge_rand_df$HasMotif <- apply(merge_rand_df[,c("Start", "End")], 1, has_motif, motif_df=fwd_motif) | apply(merge_rand_df[,c("Start", "End")], 1, has_motif, motif_df=rev_motif)

p_random <- plot_ht(merge_rand_df, energy_cutoff)
pdf(file.path(oSP, "seq_pair.random.ht.pdf"), width=30/inche_cm, height=(5+2*log2(nrow(merge_rand_df)))/inche_cm)
print(p_random)
dev.off()
cat("seq_pair.random#:", nrow(merge_rand_df))

Nanopore_enriched <- read_delim(file.path(out_dir, "Nanopore/Nanopore.enriched.tsv"), "\t", escape_double = FALSE, trim_ws = TRUE)
Nanopore_enriched <- Nanopore_enriched[, c("Start", "End")]
Nanopore_df <- inner_join(Nanopore_enriched, merge_df)
p_Nanopore_high_expr <- plot_ht(Nanopore_df, energy_cutoff)
pdf(file.path(oSP, "seq_pair.Nanopore.high_expr.ht.pdf"), width=30/inche_cm, height=(5+2*log2(nrow(Nanopore_df)))/inche_cm)
print(p_Nanopore_high_expr)
dev.off()
cat("seq_pair.Nanopore.high_expr#:", nrow(Nanopore_df))

p_Nanopore_enriched_by_gap <- plot_ht(Nanopore_df, energy_cutoff, split_group="GapType")
pdf(file.path(oSP, "seq_pair.Nanopore.high_expr.by_gap.ht.pdf"), width=30/inche_cm, height=(5+2*log2(nrow(Nanopore_df)))/inche_cm)
print(p_Nanopore_enriched_by_gap)
dev.off()

NGS_enriched <- read_delim(file.path(out_dir, "NGS/NGS.enriched.tsv"), "\t", escape_double = FALSE, trim_ws = TRUE)
NGS_enriched <- NGS_enriched[, c("Start", "End")]
NGS_df <- inner_join(NGS_enriched, merge_df)
p_NGS_high_expr <- plot_ht(NGS_df, energy_cutoff)
pdf(file.path(oSP, "seq_pair.NGS.high_expr.ht.pdf"), width=30/inche_cm, height=(5+2*log2(nrow(NGS_enriched)))/inche_cm)
print(p_NGS_high_expr)
dev.off()
cat("seq_pair.NGS.high_expr#:", nrow(NGS_df))

p_NGS_enriched_by_gap <- plot_ht(NGS_df, energy_cutoff, split_group="GapType")
pdf(file.path(oSP, "seq_pair.NGS.high_expr.by_gap.ht.pdf"), width=30/inche_cm, height=(5+2*log2(nrow(NGS_df)))/inche_cm)
print(p_NGS_enriched_by_gap)
dev.off()

all_consistent_gap <- read_delim(file.path(out_dir, "NGS_Nanopore/all_consistent_gap.tsv"), "\t", escape_double = FALSE, trim_ws = TRUE)
all_consistent_gap <- all_consistent_gap[, c("Start", "End")]
consistent_df <- inner_join(all_consistent_gap, merge_df)
p_consistent <- plot_ht(consistent_df, energy_cutoff)
pdf(file.path(oSP, "seq_pair.all_consistent.ht.pdf"), width=30/inche_cm, height=(5+2*log2(nrow(consistent_df)))/inche_cm)
print(p_consistent)
dev.off()
cat("seq_pair.seq_pair.all_consistent#:", nrow(consistent_df))


p_consistent_by_Gap <- plot_ht(consistent_df, energy_cutoff, split_group="GapType")
pdf(file.path(oSP, "seq_pair.all_consistent.by_gap.ht.pdf"), width=30/inche_cm, height=(5+2*log2(nrow(consistent_df)))/inche_cm)
print(p_consistent_by_Gap)
dev.off()

p_consistent_by_subgenome <- plot_ht(consistent_df, energy_cutoff, split_group="SubgenomeType")
pdf(file.path(oSP, "seq_pair.all_consistent.by_subgenome.ht.pdf"), width=30/inche_cm, height=(5+2*log2(nrow(consistent_df)))/inche_cm)
print(p_consistent_by_subgenome)
dev.off()

p_consistent_by_readnum <- plot_ht(consistent_df, energy_cutoff, order_indx = order(consistent_df$ReadNum, decreasing = TRUE))
pdf(file.path(oSP, "seq_pair.all_consistent.by_readnum.ht.pdf"), width=30/inche_cm, height=(5+2*log2(nrow(consistent_df)))/inche_cm)
print(p_consistent_by_readnum)
dev.off()


consistent_leader_df <- consistent_df[consistent_df$Group==consistent_df$Name, ]
p_leader_consistent <- plot_ht(consistent_leader_df, energy_cutoff)
pdf(file.path(oSP, "seq_pair.leader_consistent.ht.pdf"), width=30/inche_cm, height=(5+2*log2(nrow(consistent_leader_df)))/inche_cm)
print(p_leader_consistent)
dev.off()
cat("seq_pair.seq_pair.leader_consistent#:", nrow(consistent_leader_df))

p_leader_consistent_by_Gap <- plot_ht(consistent_leader_df, energy_cutoff, split_group="GapType")
pdf(file.path(out_dir, "seq_pair", "seq_pair.leader_consistent.by_gap.ht.pdf"), width=30/inche_cm, height=(5+2*log2(nrow(consistent_leader_df)))/inche_cm)
print(p_leader_consistent_by_Gap)
dev.off()

p_leader_consistent_by_subgenome <- plot_ht(consistent_leader_df, energy_cutoff, split_group="SubgenomeType")
pdf(file.path(out_dir, "seq_pair", "seq_pair.leader_consistent.by_subgenome.ht.pdf"), width=30/inche_cm, height=(5+2*log2(nrow(consistent_leader_df)))/inche_cm)
print(p_leader_consistent_by_subgenome)
dev.off()

p_leader_consistent_by_readnum <- plot_ht(consistent_leader_df, energy_cutoff, order_indx = order(consistent_leader_df$ReadNum, decreasing = TRUE))
pdf(file.path(out_dir, "seq_pair", "seq_pair.leader_consistent.by_readnum.ht.pdf"), width=30/inche_cm, height=(5+2*log2(nrow(consistent_leader_df)))/inche_cm)
print(p_leader_consistent_by_readnum)
dev.off()

consistent_leader_leaderSequence_df <- consistent_leader_df[consistent_leader_df$GapType=="LeaderSequence",]
p_leader_consistent_leaderSequence_by_Gap <- plot_ht(consistent_leader_leaderSequence_df, energy_cutoff)
pdf(file.path(out_dir, "seq_pair", "seq_pair.leader_consistent.LeaderSequence.ht.pdf"), width=30/inche_cm, height=(5+2*log2(nrow(consistent_leader_leaderSequence_df)))/inche_cm)
print(p_leader_consistent_leaderSequence_by_Gap)
dev.off()

## Ending pairing classfication in jumping events using +genome as template: paired or unpaired state 
group_pair_df <- merge_df %>% group_by(Group) %>% filter(n() > 1, sum(ReadNum) > 3000)
DR_AR_pair_DR_mat <- pair2mat(group_pair_df$DR_AR_pair_DR, group_pair_df$Name, group_pair_df$has_DR_AR_pair, group_pair_df$DR_AR_pair_Energy, energy_cutoff)
DR_AR_pair_AR_mat <- pair2mat(group_pair_df$DR_AR_pair_AR, group_pair_df$Name, group_pair_df$has_DR_AR_pair, group_pair_df$DR_AR_pair_Energy, energy_cutoff)
group_pair_df$EndPair <- DR_AR_pair_DR_mat[, 1] & DR_AR_pair_AR_mat[, ncol(DR_AR_pair_AR_mat)] # if both in paired state

group_pair_df <- group_pair_df %>% group_by(Group) %>% filter(sum(EndPair)>=5, sum(!EndPair)>=5)
write_tsv(group_pair_df, file.path(out_dir, "seq_pair", "EndPair.SourceData.tsv"))
group_pair_df$EndPair <- factor(group_pair_df$EndPair, levels = c(TRUE, FALSE), labels = c("Paired", "Unpaired"))

group_pair_leader_df <- group_pair_df[group_pair_df$Name == group_pair_df$Group, ]
group_pair_leader_df$PairType <- "DL-AL"
group_pair_leader_df$PairType[(group_pair_leader_df$DL_AL_pair_Energy - group_pair_leader_df$DR_AR_pair_Energy) > 0] <- "DR-AR"
group_pair_leader_df$Label <- sprintf(
  "Gap: %s\nGap type: %s\nSubgenome type: %s\nPaired type: %s\nEvent number: %d\nRead number: %d", 
  group_pair_leader_df$Group, 
  group_pair_leader_df$GapType, 
  group_pair_leader_df$SubgenomeType, 
  group_pair_leader_df$PairType,
  group_pair_leader_df$GroupSize,
  group_pair_leader_df$GroupReadNum
  )
group_pair_leader_df <- group_pair_leader_df[order(group_pair_leader_df$GroupReadNum), ]
group_pair_leader_df <- group_pair_leader_df[, c("Group", "Label")]
group_pair_leader_df$Label <- factor(group_pair_leader_df$Label, levels = rev(group_pair_leader_df$Label))

plot_df <- group_pair_df[, c("Group", "EndPair", "ReadNum")]
plot_df <- left_join(plot_df, group_pair_leader_df)
p <- ggplot(plot_df, aes(x=EndPair, y=log10(ReadNum))) +
  geom_boxplot(aes(fill=EndPair), size=0.2, outlier.color = NA) +
  geom_point(color="black", position = position_jitter(width = 0.2), size=0.5) +
  facet_wrap(~Label, nrow=1) +
  theme_bw() +
  scale_y_continuous(breaks = c(0, 3, 6), labels = c("1E0", "1E3", "1E6")) +
  labs(y="#Read") +
  theme(text = element_text(family="ArialMT", size=6),
        title = element_text(family="ArialMT", size=6),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        # legend.position = c(0.1, 0.7),
        legend.position = "none",
        panel.grid = element_blank()
  )
ggsave(file.path(out_dir, "seq_pair", "EndPair.pdf"), p, width = 15, height = 5, units = "cm")


## Motif
consistent_df$HasMotif <- factor(consistent_df$HasMotif, levels = c(TRUE, FALSE), labels = c("Include motif", "Exclude motif"))
p <- ggplot(consistent_df, aes(x=GapType, fill=HasMotif)) +
  geom_bar() +
  scale_fill_brewer(palette = "Dark2") +
  theme_bw() +
  labs(y="#Junctions") +
  theme(text = element_text(family="ArialMT", size=6),
        title = element_text(family="ArialMT", size=6),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(color = "black"),
        legend.key.size = unit(4, "mm"),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(family="ArialMT", size=6),
        panel.grid = element_blank()
  )
ggsave(file.path(out_dir, "seq_pair", "Motif.cnt.by_gap.pdf"), p, width = 5, height = 5, units = "cm")

consistent_leader_df <- consistent_df[consistent_df$GapType=="LeaderSequence",]
p <- ggplot(consistent_leader_df, aes(x=SubgenomeType, fill=HasMotif)) +
  geom_bar() +
  scale_fill_brewer(palette = "Dark2") +
  theme_bw() +
  labs(y="#Junctions") +
  theme(text = element_text(family="ArialMT", size=6),
        title = element_text(family="ArialMT", size=6),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(color = "black"),
        legend.key.size = unit(4, "mm"),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(family="ArialMT", size=6),
        panel.grid = element_blank()
  )
ggsave(file.path(out_dir, "seq_pair", "Motif.cnt.by_gap.consistent.pdf"), p, width = 5.5, height = 4, units = "cm")

merge_df$IsConsistent <- merge_df$Name %in% consistent_df$Name
merge_df$PairType <- NA
merge_df$PairType[merge_df$DL_AL_pair_Energy < merge_df$DR_AR_pair_Energy] <- "DL-AL"
merge_df$PairType[merge_df$DL_AL_pair_Energy >= merge_df$DR_AR_pair_Energy] <- "DR-AR"
merge_df$PairType[abs(merge_df$DL_AL_pair_Energy >= merge_df$DR_AR_pair_Energy)] <- "NC"
all(merge_df$has_DL_AL_pair)
all(merge_df$has_DR_AR_pair)
write_tsv(merge_df, file.path(out_dir, "seq_pair", "Pair.info.all.tsv"))
