library(readr)
library(dplyr)
library(ggplot2)
library(BiocParallel)
library(ComplexHeatmap)
library(RColorBrewer)
library(reshape2)
library(igraph)
library(circlize)

data_dir <- file.path("data", "other_data")
out_dir <- file.path("analysis", "subgenome_other_data")
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
  
  cb_palette <- c("#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e",
                  "#4aef7b", "#e86502", "#9ed84e", "#39ba30", "#6ad157", "#8249aa", "#99db27", "#e07233", "#ff523f",
                  "#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
                  "#1a918f", "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
                  "#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
                  "#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5","#925bea", "#63ff4f")
  
  SubgenomeType_col <- c(cb_palette[1:(length(levels(merged_df_row_clu$SubgenomeType))-1)], "grey70")
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



for(s in c("SARS", "MERS")){
  js_data_dir <- file.path(data_dir, "virus_JS", s)
  NGS_JS_rep1 <- load_NGS_js(file.path(js_data_dir, sprintf("%s.rep1_jd.csv", s)))
  NGS_JS_rep2 <- load_NGS_js(file.path(js_data_dir, sprintf("%s.rep2_jd.csv", s)))
  both_NGS_JS <- inner_join(NGS_JS_rep1, NGS_JS_rep2)
  
  l_df <- load_pair(file.path(data_dir, "virus_pairing", s, sprintf("%s.rep1.lft.flank.csv", s)), 
                    file.path(data_dir, "virus_pairing", s, sprintf("%s.rep1.lft_hybcons.csv", s)), both_NGS_JS)
  r_df <- load_pair(file.path(data_dir, "virus_pairing", s, sprintf("%s.rep1.rgt.flank.csv", s)), 
                    file.path(data_dir, "virus_pairing", s, sprintf("%s.rep1.rgt_hybcons.csv", s)), both_NGS_JS)
  
  df <- inner_join(l_df[, c("Start", "End")], r_df[, c("Start", "End")])
  l_df <- left_join(df, l_df)
  r_df <- left_join(df, r_df)
  names(l_df)[4:9] <- c("DL_seq", "AL_seq", "DL_AL_pair_DL", "DL_AL_pair_AL", "DL_AL_pair_Energy", "has_DL_AL_pair")
  names(r_df)[4:9] <- c("DR_seq", "AR_seq", "DR_AR_pair_DR", "DR_AR_pair_AR", "DR_AR_pair_Energy", "has_DR_AR_pair")
  
  LEADER_BOUNDARY <- 100
  ORF1ab_BOUNDARY <- 20000
  
  group_df <- creat_group(l_df, resolution=5)
  group_df$GapType <- "Other"
  group_df$GapType[group_df$Start <= LEADER_BOUNDARY] <- "LeaderSequence"
  group_df$GapType[(group_df$Start > LEADER_BOUNDARY) & (group_df$Start <= ORF1ab_BOUNDARY) & (group_df$End > ORF1ab_BOUNDARY)] <- "ORF1ab_S2N"
  group_df$GapType[(group_df$Start > ORF1ab_BOUNDARY) & (group_df$End > ORF1ab_BOUNDARY)] <- "S2N_S2N"
  group_df$GapType <- factor(group_df$GapType, levels = c("LeaderSequence", "ORF1ab_S2N", "S2N_S2N", "Other"))
  
  gene_ref_df <- read_delim(file.path(data_dir, sprintf("%s.ref.tsv", s)),
                         "\t", escape_double = FALSE, trim_ws = TRUE)
  gene_ref_df$Mid <- (gene_ref_df$Start + gene_ref_df$End) / 2
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
  
  l_rep2_df <- load_pair(file.path(data_dir, "virus_pairing", s, sprintf("%s.rep2.lft.flank.csv", s)), 
                    file.path(data_dir, "virus_pairing", s, sprintf("%s.rep2.lft_hybcons.csv", s)), both_NGS_JS)
  l_rep2_df <- left_join(group_df, group_df)
  merge_df$ReadNum <- merge_df$ReadNum + l_rep2_df$ReadNum
  
  leader_js_df <- group_df[group_df$GapType=="LeaderSequence" & group_df$End>ORF1ab_BOUNDARY,]
  end_cnt <- leader_js_df %>% group_by(SubgenomeType, LeaderEnd) %>% 
    summarise(ReadNum=sum(ReadNum)) %>% 
    group_by(SubgenomeType) %>% 
    arrange(-ReadNum, .by_group = TRUE) %>%
    mutate(Major=ReadNum==max(ReadNum))
  y_indx <- data.frame(SubgenomeType=c(gene_ref_df$Name[gene_ref_df$Name!="orf1ab"], "Other"))
  y_indx$y_indx <- 1:nrow(y_indx)
  end_cnt <- left_join(end_cnt, y_indx)
  end_cnt$y_end <- end_cnt$y_indx - 0.3
  end_cnt$y_end[end_cnt$Major] <- end_cnt$y_indx[end_cnt$Major] - 0.5
  
  y_indx$Name <- y_indx$SubgenomeType
  tmp_gene_ref_df <- gene_ref_df[gene_ref_df$Name!="orf1ab",] 
  tmp_gene_ref_df <- left_join(tmp_gene_ref_df, y_indx)
  
  min_site_df <- end_cnt %>% group_by(LeaderEnd, SubgenomeType, y_indx) %>% summarise(min_site=min(LeaderEnd))
  
  ## AMGAAC motif
  TRScore <- read_delim(file.path(data_dir, sprintf("%s.TRS.bed", s)), "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  names(TRScore) <- c("Chrom", "Start", "End", "Name", "Score", 'Strand')
  TRScore <- TRScore[,c("Start", "End", "Name", "Strand")]
  TRScore <- TRScore[order(TRScore$Start),]
  
  start_num <- merge_df %>% group_by(Start) %>% summarise(Number=sum(ReadNum))
  names(start_num)[1] <- "Position"
  end_num <- merge_df %>% group_by(End) %>% summarise(Number=sum(ReadNum))
  names(end_num)[1] <- "Position"
  
  
  
  fwd_motif <- TRScore[TRScore$Strand=="+",]
  rev_motif <- TRScore[TRScore$Strand=="-",]
  
  has_motif <- function(line, motif_df, expand=20){
    donor_pos <- line[1]
    acceptor_pos <- line[2]
    donor_match <- (motif_df$Start >= (donor_pos-expand)) & (motif_df$End <= (donor_pos+expand))
    acceptor_match <- (motif_df$Start >= (acceptor_pos-expand)) & (motif_df$End <= (acceptor_pos+expand))
    pair_match <- any(donor_match) & any(acceptor_match)
    return(any(pair_match))
  }
  merge_df$HasMotif <- apply(merge_df[,c("Start", "End")], 1, has_motif, motif_df=fwd_motif) | apply(merge_df[,c("Start", "End")], 1, has_motif, motif_df=rev_motif)
  write_tsv(merge_df, file.path(oSP, sprintf("seq_pair.SourceData.%s.tsv", s)))
  
  
  
  # Analyze MFE
  merge_df$Name <- sprintf("%d-%d", merge_df$Start, merge_df$End)
  energy_cutoff <- 0

  top_js <- merge_df %>% group_by(SubgenomeType) %>% filter(ReadNum==max(ReadNum))
  top_js <- top_js[!top_js$SubgenomeType %in% c("Other", "orf1ab"),]
  top_js$MFE <- apply(top_js[,c("DL_AL_pair_Energy", "DR_AR_pair_Energy")], 1, min)
  top_js <- top_js[order(top_js$End),]
  top_js$SubgenomeType <- factor(top_js$SubgenomeType, levels = top_js$SubgenomeType)

  rank_df <- merge_df
  rank_df$MFE <- apply(merge_df[,c("DL_AL_pair_Energy", "DR_AR_pair_Energy")], 1, min)
  rank_df <- rank_df %>% group_by(SubgenomeType) %>% mutate(IsTop=ReadNum==max(ReadNum))
  rank_df <- rank_df[rank_df$SubgenomeType%in%gene_ref_df$Name[2:nrow(gene_ref_df)],]
  rank_df <- rank_df[order(rank_df$ReadNum,decreasing = TRUE),]
  rank_df$Order <- 1:nrow(rank_df)
  rank_df$SubgenomeType <- factor(rank_df$SubgenomeType, levels = top_js$SubgenomeType)
  
  ################################################
  ## Figure S7C
  
  rank_df$Class <- "<100"
  rank_df$Class[rank_df$ReadNum>=100] <- "100~1000"
  rank_df$Class[rank_df$ReadNum>1000] <- ">1000"
  rank_df$Class <- factor(rank_df$Class, levels = c("<100", "100~1000", ">1000"))
  rank_df_info <- rank_df %>% group_by(Class) %>% summarise(n=n())
  rank_df_info$Label <- sprintf("n=%d", rank_df_info$n)
  
  p <- ggplot(rank_df, aes(x=Class, y=MFE, fill=Class)) +
    geom_boxplot(outlier.color = NA, size=0.2) +
    geom_text(data=rank_df_info, aes(y=-5, label=Label), size=1.8) +
    scale_fill_brewer(palette = "Dark2") +
    labs(x="#Nanopore reads", y="MFE (kcal/mol)") +
    theme_bw() +
    theme(text = element_text(family="ArialMT", size=5),
          title = element_text(family="ArialMT", size=5),
          axis.text = element_text(color = "black"),
          legend.position = "none",
          panel.grid = element_blank()
    )
  ggsave(file.path(oSP, sprintf("MFE.ReadNumMFE.boxplot.%s.pdf", s)), p, width = 3, height = 4, units = "cm")
  
  ################################################
  ## Figure S7B
  
  inche_cm <- 2.54
  p_merge_by_gap <- plot_ht(merge_df[merge_df$ReadNum>=10,], energy_cutoff, split_group="GapType")
  pdf(file.path(oSP, sprintf("seq_pair.all.by_gap.ht.%s.pdf", s)), width=30/inche_cm, height=(5+2*log2(nrow(merge_df)))/inche_cm)
  print(p_merge_by_gap)
  dev.off()
  
  
  leader_df <- merge_df[merge_df$SubgenomeType %in% gene_ref_df$Name[2:nrow(gene_ref_df)],]
  leader_df$R_strong <- leader_df$DR_AR_pair_Energy < leader_df$DL_AL_pair_Energy
  leader_DR_AR_pair_DR_mat <- pair2mat(leader_df$DR_AR_pair_DR, leader_df$Name, leader_df$has_DR_AR_pair, leader_df$DR_AR_pair_Energy, energy_cutoff)
  leader_DR_AR_pair_AR_mat <- pair2mat(leader_df$DR_AR_pair_AR, leader_df$Name, leader_df$has_DR_AR_pair, leader_df$DR_AR_pair_Energy, energy_cutoff)
  leader_DL_AL_pair_DL_mat <- pair2mat(leader_df$DL_AL_pair_DL, leader_df$Name, leader_df$has_DL_AL_pair, leader_df$DL_AL_pair_Energy, energy_cutoff)
  leader_DL_AL_pair_AL_mat <- pair2mat(leader_df$DL_AL_pair_AL, leader_df$Name, leader_df$has_DL_AL_pair, leader_df$DL_AL_pair_Energy, energy_cutoff)
  R_EndPair <- leader_DR_AR_pair_DR_mat[, 1] & leader_DR_AR_pair_AR_mat[, ncol(leader_DR_AR_pair_AR_mat)] & leader_df$R_strong
  L_EndPair <- leader_DL_AL_pair_DL_mat[, ncol(leader_DL_AL_pair_DL_mat)] & leader_DL_AL_pair_AL_mat[, 1] & (!leader_df$R_strong)
  leader_df$EndPair <- R_EndPair | L_EndPair # if both in paired state
  leader_df$MFE <- apply(leader_df[,c("DL_AL_pair_Energy", "DR_AR_pair_Energy")], 1, min)
  
  leader_df <- leader_df %>% group_by(Group, R_strong) %>% filter(sum(EndPair)>=1, sum(!EndPair)>=1, max(ReadNum)>100)
  leader_df <- leader_df[leader_df$R_strong,]
  
  find_endpair_js <- function(block, max_unmatch=2){
    end_pair_df <- block[block$EndPair,]
    no_end_pair_df <- block[!block$EndPair,]
    end_pair_name_li <- c()
    no_end_pair_name_li <- c()
    end_pair_DR_AR_pair_DR_mat <- pair2mat(end_pair_df$DR_AR_pair_DR, end_pair_df$Name, end_pair_df$has_DR_AR_pair, end_pair_df$DR_AR_pair_Energy, energy_cutoff)
    end_pair_DR_AR_pair_AR_mat <- pair2mat(end_pair_df$DR_AR_pair_AR, end_pair_df$Name, end_pair_df$has_DR_AR_pair, end_pair_df$DR_AR_pair_Energy, energy_cutoff)
    no_end_pair_DR_AR_pair_DR_mat <- pair2mat(no_end_pair_df$DR_AR_pair_DR, no_end_pair_df$Name, no_end_pair_df$has_DR_AR_pair, no_end_pair_df$DR_AR_pair_Energy, energy_cutoff)
    no_end_pair_DR_AR_pair_AR_mat <- pair2mat(no_end_pair_df$DR_AR_pair_AR, no_end_pair_df$Name, no_end_pair_df$has_DR_AR_pair, no_end_pair_df$DR_AR_pair_Energy, energy_cutoff)
    
    is_match <- function(indx1, indx2, max_unmatch=2){
      max_pair_num <- max(length(indx1), length(indx2))
      for(shift in -10:10){
        tmp_indx <- indx1 + shift
        overlap_num <- length(intersect(tmp_indx, indx2))
        if((max_pair_num - overlap_num) <= max_unmatch){
          return(TRUE)
        }
      }
      return(FALSE)
    }
    for(i in 1:nrow(end_pair_df)){
      for(j in 1:nrow(no_end_pair_df)){
        end_pair_DR_indx <- which(end_pair_DR_AR_pair_DR_mat[i,])
        no_end_pair_DR_indx <- which(no_end_pair_DR_AR_pair_DR_mat[j,])
        end_pair_AR_indx <- which(end_pair_DR_AR_pair_AR_mat[i,])
        no_end_pair_AR_indx <- which(no_end_pair_DR_AR_pair_AR_mat[j,])
        DR_match <- is_match(end_pair_DR_indx, no_end_pair_DR_indx, max_unmatch)
        AR_match <- is_match(end_pair_AR_indx, no_end_pair_AR_indx, max_unmatch)
        if(DR_match & AR_match){
          end_pair_name_li <- c(end_pair_name_li, end_pair_df$Name[i])
          no_end_pair_name_li <- c(no_end_pair_name_li, no_end_pair_df$Name[j])
        }
      }
    }
    if(length(end_pair_name_li)==0){
      return(data.frame())
    }
    edge_df <- data.frame(from=end_pair_name_li, to=no_end_pair_name_li)
    edge_df <- edge_df[order(edge_df$from, edge_df$to), ]
    subgenome_graph <- graph_from_edgelist(as.matrix(edge_df), directed = FALSE)
    graph_cluster <- components(subgenome_graph, "strong")
    
    res <- data.frame(
      Name=names(graph_cluster$membership), 
      PairClu=graph_cluster$membership
    )
    res <- left_join(res, block)
    return(res)
  }
  
  pair_js_df <- plyr::ddply(leader_df, "Group", find_endpair_js)
  pair_js_subgroup_info <- pair_js_df %>% 
    group_by(Group, PairClu) %>% 
    summarise(
      SubGroup=sprintf("(%d~%d)-(%d~%d)", min(Start), max(Start), min(End), max(End)), 
      IsMatch=EndPair[which.max(ReadNum)]
    )
  pair_js_df <- left_join(pair_js_df, pair_js_subgroup_info)
  
  write_tsv(pair_js_df, file.path(out_dir, "seq_pair", sprintf("EndPair.sourceData.%s.tsv", s)))
  pair_js_df$EndPair <- factor(pair_js_df$EndPair, levels = c(TRUE, FALSE), labels = c("Paired", "Unpaired"))
  
  ################################################
  ## Figure S7D
  p <- plot_ht(pair_js_df, energy_cutoff, order_indx = order(pair_js_df$SubGroup, -pair_js_df$ReadNum), split_group = "SubGroup")
  pdf(file.path(out_dir, "seq_pair", sprintf("EndPair.ht.%s.pdf", s)), width=30/inche_cm, height=10/inche_cm)
  print(p)
  dev.off()
}


