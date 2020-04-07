library(readr)
library(dplyr)
library(igraph)
library(ggplot2)
library(reshape2)
library(circlize)
library(RColorBrewer)
library(rtracklayer)


data_dir <- file.path("data", "subgenome")
out_dir <- file.path("analysis", "subgenome")
if (!dir.exists(out_dir)) dir.create(out_dir)

fl_cutoff <- 45
chrom_size <- 29891

gene_ref <- read_delim(file.path(data_dir, "genome/WIV04.bed"), "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
gene_ref_df <- data.frame(Name=gene_ref$X4, Start=gene_ref$X2, End=gene_ref$X3)
gene_ref_df$Mid <- (gene_ref_df$Start + gene_ref_df$End) / 2
gene_ref_df <- gene_ref_df[order(gene_ref_df$Start), ]

all_consistent_gap <- read_delim(file.path(out_dir, "NGS_Nanopore", "all_consistent_gap.tsv"), 
                                 "\t", escape_double = FALSE, trim_ws = TRUE)
nanopore_bed_rep1 <- read_delim(file.path(data_dir, "nanopore", "data", "WIV04.rep1.bed"), 
                                "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
nanopore_bed_rep2 <- read_delim(file.path(data_dir, "nanopore", "data", "WIV04.rep2.bed"), 
                                "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
# Keep FL reads
nanopore_rep1_fl <- nanopore_bed_rep1$X4[(nanopore_bed_rep1$X2 < fl_cutoff) & (nanopore_bed_rep1$X3 > (chrom_size-fl_cutoff))]
nanopore_rep2_fl <- nanopore_bed_rep2$X4[(nanopore_bed_rep2$X2 < fl_cutoff) & (nanopore_bed_rep2$X3 > (chrom_size-fl_cutoff))]

nanopore_rep1_fl_genome_read <- nanopore_bed_rep1[nanopore_bed_rep1$X10==1 & (nanopore_bed_rep1$X4 %in% nanopore_rep1_fl),]
nanopore_rep2_fl_genome_read <- nanopore_bed_rep2[nanopore_bed_rep2$X10==1 & (nanopore_bed_rep2$X4 %in% nanopore_rep2_fl),]
nanopore_rep1_fl_genome_read_num <- nrow(nanopore_rep1_fl_genome_read)
nanopore_rep2_fl_genome_read_num <- nrow(nanopore_rep2_fl_genome_read)
nanopore_fl_genome_read_num <- nanopore_rep1_fl_genome_read_num + nanopore_rep2_fl_genome_read_num

# Filter out reads with multiple mapping hits
nanopore_multiple_rep1 <- nanopore_bed_rep1 %>% group_by(X4) %>% filter(n()>1) %>% summarise(MapNum=n())
nanopore_multiple_rep2 <- nanopore_bed_rep2 %>% group_by(X4) %>% filter(n()>1) %>% summarise(MapNum=n())
nanopore_rep1_fl <- setdiff(nanopore_rep1_fl, nanopore_multiple_rep1$X4)
nanopore_rep2_fl <- setdiff(nanopore_rep2_fl, nanopore_multiple_rep2$X4)

## Only JS within FL reads will be used
nanopore_JS_rep1 <- read_delim(file.path(data_dir, "nanopore", "data", "WIV04.rep1.gap.bed"), 
                               "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
names(nanopore_JS_rep1) <- c("Chrom", "Start", "End", "Name", "Score", "Strand")
nanopore_JS_rep1 <- nanopore_JS_rep1[nanopore_JS_rep1$Name %in% nanopore_rep1_fl, ]
nanopore_JS_rep1$Source <- "Rep1"

nanopore_JS_rep2 <- read_delim(file.path(data_dir, "nanopore", "data", "WIV04.rep2.gap.bed"), 
                               "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
names(nanopore_JS_rep2) <- c("Chrom", "Start", "End", "Name", "Score", "Strand")
nanopore_JS_rep2 <- nanopore_JS_rep2[nanopore_JS_rep2$Name %in% nanopore_rep2_fl, ]
nanopore_JS_rep2$Source <- "Rep2"

merged_nanopore_JS <- rbind(nanopore_JS_rep1, nanopore_JS_rep2)
merged_nanopore_JS <- merged_nanopore_JS[order(merged_nanopore_JS$Start, merged_nanopore_JS$End), ]
merged_nanopore_JS <- merged_nanopore_JS[!duplicated(merged_nanopore_JS), ]
all_consistent_JuncList <- sprintf("%s-%s", all_consistent_gap$Start, all_consistent_gap$End)

subgenome_df <- merged_nanopore_JS %>% 
  group_by(Name, Source) %>% 
  filter(all(sprintf("%s-%s", Start, End) %in% all_consistent_JuncList)) %>% 
  summarise(JuncList=paste(sprintf("%s-%s", Start, End), collapse = " "), JuncNum = n()) %>% 
  group_by(JuncList, JuncNum) %>% 
  summarise(Number = n())
subgenome_df <- subgenome_df[order(subgenome_df$Number, decreasing = T), ]
head(subgenome_df)
cat("Number of full-length subgenomes: ", nrow(subgenome_df))

# Save full-length reads supported subgenomes
oSG <- file.path(out_dir, "subgenome")
if (!dir.exists(oSG)) dir.create(oSG)
write_tsv(subgenome_df, file.path(oSG, "subgenome.raw.tsv"))

#' Extract junctions in subgenomes
juncStr2df <- function(x) {
  df <- as.data.frame(t(as.data.frame(strsplit(x, "-"))))
  names(df) <- c("Start", "End")
  df$Start <- as.integer(as.character(df$Start))
  df$End <- as.integer(as.character(df$End))
  df <- df[order(df$Start), ]
  rownames(df) <- NULL
  return(df)
}

#' Cluster similar subgenomes with the same jumping structure
group_subgenome <- function(subgenome_df, resolution=5) {
  junc_df_li <- lapply(strsplit(subgenome_df$JuncList, " "), juncStr2df)
  
  #' require same #junctions and position differces of all junctions less than cutoff
  is_one_group <- function(junc_df1, junc_df2, resolution=5) {
    if(nrow(junc_df1) != nrow(junc_df2)) {
      return(FALSE)
    } 
    diff_junc <- junc_df1 - junc_df2
    return(all(abs(diff_junc) <= resolution))
  }
  indx1 <- c()
  indx2 <- c()
  for(i in 1:nrow(subgenome_df)) {
    for(j in 1:i) {
      can_merge <- is_one_group(junc_df_li[[i]], junc_df_li[[j]], resolution)
      if(can_merge){
        indx1 <- c(indx1, i)
        indx2 <- c(indx2, j)
      }
    }
  }
  
  # create connecting graph
  edge_df <- data.frame(from=indx1, to=indx2)
  edge_df <- edge_df[order(edge_df$from, edge_df$to), ]
  subgenome_graph <- graph_from_edgelist(as.matrix(edge_df), directed = FALSE)
  graph_clustar <- components(subgenome_graph, "strong")
  
  # merge subgenomes by group
  LEADER_BOUNDARY <- 100
  ORF1ab_BOUNDARY <- 20000
  merge_subgenome <- function(indx_li, subgenome_df, junc_df_li) {
    selected_df <- subgenome_df[indx_li, c("JuncList", "JuncNum", "Number")]
    selected_df <- selected_df[order(selected_df$Number, decreasing = TRUE), ]
    junc_df <- junc_df_li[[indx_li[1]]]
    
    # locate the first downstream complete gene of the first junction
    idx <- which(min(junc_df$End) < gene_ref_df$Start)
    if (length(idx) == 0) {
      sgtype <- "Other"
    } else {
      sgtype <- gene_ref_df$Name[min(idx)] 
    }
    
    res <- data.frame(
      LeaderJuncList = selected_df$JuncList[1],  # core subgenome within one group by number of reads
      JuncNum = selected_df$JuncNum[1], 
      EventNumber = length(indx_li), # number of similar subgenomes in one group
      LeaderReadNumber = selected_df$Number[1], 
      GroupReadNumber = sum(selected_df$Number),
      LeaderSequence = sum(junc_df$Start <= LEADER_BOUNDARY), # group starting in leader sequence
      ORF1ab_S2N = sum(junc_df$Start > LEADER_BOUNDARY & junc_df$Start <= ORF1ab_BOUNDARY & junc_df$End > ORF1ab_BOUNDARY),
      S2N_S2N = sum(junc_df$Start > ORF1ab_BOUNDARY),
      Type = sgtype
      )
    return(res)
  }
  merged_group <- do.call(rbind, lapply(
    groups(graph_clustar), merge_subgenome, subgenome_df = subgenome_df, junc_df_li = junc_df_li))
  merged_group <- merged_group[order(merged_group$GroupReadNumber, decreasing = TRUE), ]
  return(merged_group)
}

subgenome_group_df <- group_subgenome(subgenome_df, resolution = 5)
subgenome_group_df$Type <- as.character(subgenome_group_df$Type)
head(subgenome_group_df)
write_tsv(subgenome_group_df, file.path(oSG, "subgenome.group.tsv"))

# Summarize the subgenome structure patterns
junc_class_info <- subgenome_group_df %>% 
  group_by(JuncNum, LeaderSequence, ORF1ab_S2N, S2N_S2N) %>%
  summarise(GroupNumber=n(), EventNumber=sum(EventNumber), ReadNumber=sum(GroupReadNumber))
junc_class_info

junc_class_info <- junc_class_info[order(junc_class_info$ReadNumber, decreasing = TRUE), ]
junc_class_info$Structure <- paste(junc_class_info$LeaderSequence, junc_class_info$ORF1ab_S2N, junc_class_info$S2N_S2N)
junc_class_info$Structure <- factor(
  junc_class_info$Structure, 
  levels = rev(c(
    "1 0 0", 
    "1 0 1", 
    "0 1 0", 
    "1 0 2",
    "0 0 1")), 
  labels = rev(c(
    "#----------------#####",
    "#----------------#--##",
    "####-------------#####",
    "#----------------#-#-#",
    "##################--##"))
)
write_tsv(junc_class_info, file.path(out_dir, "subgenome", "subgenome.group.stat.tsv"))

# Plot the pattern and numbers
junc_class_info <- junc_class_info[, c("Structure", "GroupNumber", "EventNumber", "ReadNumber")]
junc_class_df <- melt(junc_class_info, "Structure", variable.name = "Type", value.name = "Number")
junc_class_df$Type <- factor(junc_class_df$Type, 
                             levels = c("GroupNumber", "EventNumber", "ReadNumber"), 
                             labels = c("#Cluster", "#Event", "#Read"))
p <- ggplot(junc_class_df, aes(x=Structure, y=Number, fill=Structure)) +
  geom_bar(stat="identity", position="dodge", width = 0.7) +
  geom_text(aes(label=Number), size=1.8, hjust=1.1) +
  theme_bw() +
  scale_fill_brewer(palette="Set1") +
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
ggsave(file.path(oSG, "subgenome.class.pdf"), p, width = 8, height = 2.5, units = "cm")


# Summarize the first complete gene downstream of the jumping
gene_type_info <- subgenome_group_df %>% 
  group_by(Type) %>%
  summarise(GroupNumber=n(), EventNumber=sum(EventNumber), ReadNumber=sum(GroupReadNumber))
gene_type_info$Type <- factor(gene_type_info$Type, levels = c("Other", as.character(gene_ref_df$Name)))
gene_type_info <- gene_type_info[order(gene_type_info$Type, decreasing = TRUE),]
gene_type_info$CumReadNum <- NA
gene_type_info$CumReadNum[1:9] <- rev(cumsum(rev(gene_type_info$ReadNumber[1:9])))
gene_type_info
write_tsv(gene_type_info, file.path(oSG, "subgenome.type.stat.tsv"))


# Plot the groups and numbers
junc_class_df <- melt(gene_type_info, "Type", variable.name = "Stat", value.name = "Number")
junc_class_df$Stat <- factor(junc_class_df$Stat, 
                             levels = c("GroupNumber", "EventNumber", "ReadNumber", "CumReadNum"), 
                             labels = c("#Cluster", "#Event", "#Read", "Cumulative #Read"))

SubgenomeType_col <- rev(c(brewer.pal(length(levels(junc_class_df$Type))-1, "RdBu"), "grey70"))
names(SubgenomeType_col) <- levels(junc_class_df$Type)

p <- ggplot(junc_class_df, aes(x=Type, y=Number, fill=Type)) +
  geom_bar(stat="identity", position="dodge", width = 0.7, size=0.3, color="black") +
  geom_text(aes(label=Number), size=1.8, hjust=1.1) +
  theme_bw() +
  scale_fill_manual(values = SubgenomeType_col) +
  facet_wrap(~Stat, scale="free_x", nrow=1) +
  coord_flip() +
  # scale_y_log10() +
  theme(text = element_text(family="ArialMT", size=6),
        title = element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        panel.grid = element_blank()
  )
ggsave(file.path(oSG, "subgenome.type.pdf"), p, width = 12, height = 4, units = "cm")

# Build data for presenting subgenome structure and statistics
build_str_df <- function(subgenome_group_df, chrom_size=29891){
  junc_df_li <- lapply(strsplit(as.character(subgenome_group_df$LeaderJuncList), " "), juncStr2df)
  exon_df <- data.frame()
  intron_df <- data.frame()
  for (i in 1:nrow(subgenome_group_df)) {
    junc_df <- junc_df_li[[i]]
    junc_num <- nrow(junc_df)
    junc_df$y <- i
    junc_df$Str <- subgenome_group_df$Structure[i]
    intron_df <- rbind(intron_df, junc_df)
    exon_start <- c(0, junc_df$End[junc_num])
    exon_end <- c(junc_df$Start[1], chrom_size)
    if (junc_num >1) {
      for (j in 2:junc_num) {
        exon_start <- c(exon_start, junc_df$End[j-1])
        exon_end <- c(exon_end, junc_df$Start[j])
      }
    }
    exon_df <- rbind(exon_df, data.frame(
      y=i, Start=exon_start, End=exon_end, Str=subgenome_group_df$Structure[i]))
  }
  number_df <- data.frame(
    y = 1:nrow(subgenome_group_df), 
    ReadNum = subgenome_group_df$GroupReadNumber, 
    LogReadNum = log10(subgenome_group_df$GroupReadNumber+1),
    Str = subgenome_group_df$Structure,
    Type = subgenome_group_df$Type,
    Label = sprintf("%s (%d)", subgenome_group_df$Type, subgenome_group_df$GroupReadNumber)
    )
  return(list(
    "intron" = intron_df,
    "exon" = exon_df,
    "number" = number_df
  ))
}

subgenome_group_df$Structure <- paste(
  subgenome_group_df$LeaderSequence, subgenome_group_df$ORF1ab_S2N, subgenome_group_df$S2N_S2N)
subgenome_group_df$Structure <- factor(
  subgenome_group_df$Structure, 
  levels = rev(c(
    "1 0 0", 
    "1 0 1", 
    "0 1 0", 
    "1 0 2",
    "0 0 1")), 
  labels = rev(c(
    "#----------------#####",
    "#----------------#--##",
    "####-------------#####",
    "#----------------#-#-#",
    "##################--##"))
)
subgenome_group_df <- subgenome_group_df[
  order(subgenome_group_df$GroupReadNumber, subgenome_group_df$Structure, decreasing = TRUE), ]
num_start <- chrom_size * 1.1

# Plot
plot_str_li <- build_str_df(subgenome_group_df, chrom_size=chrom_size)
p <- ggplot() +
  geom_rect(data=gene_ref_df, aes(xmin=Start, xmax=End, ymin=-0.3, ymax=0.3), fill="black") + # gene track
  geom_text(data=gene_ref_df, aes(x=Mid, y=1, label=Name), color="black", size=1.2) +
  geom_rect(data=plot_str_li[["exon"]], aes(xmin=Start, xmax=End, ymin=-y-0.3, ymax=-y+0.3, fill=Str)) +
  geom_linerange(data=plot_str_li[["intron"]], aes(xmin=Start, xmax=End, y=-y, color=Str), size=0.3) +
  geom_rect(data=plot_str_li[["number"]], aes(xmin=num_start, xmax=num_start+1000*LogReadNum, ymin=-y-0.3, ymax=-y+0.3, fill=Str)) +
  geom_text(data=plot_str_li[["number"]], aes(x=num_start+1000*LogReadNum, y=-y, label=Label), color="black", size=1.2, hjust=-0.1) +
  geom_vline(xintercept = num_start, color="black", size=0.3) +
  scale_x_continuous(
    breaks = c(0, 15000, chrom_size, num_start+1000*log10(1+1), num_start+1000*log10(1+1000), num_start+1000*log10(1+1000000)), 
    labels = c("0", "15000", as.character(chrom_size), "1e0", "1e3", "1e6")
    ) +
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") +
  theme_bw() +
  theme(
    text = element_text(family="ArialMT", size=6),
    title = element_blank(),
    axis.text = element_text(color = "black"),
    axis.ticks.y = element_blank(),
    legend.key.size = unit(4, "mm"),
    legend.title = element_blank(),
    legend.text = element_text(family="ArialMT", size=6),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_line(color="black", size=0.3),
    axis.line.x = element_line(color="black", size=0.3),
    axis.text.y = element_blank(),
    axis.text.x = element_text(hjust=1)
  )
ggsave(file.path(oSG, "subgenome.str.pdf"), p, width = 12, height = 20, units = "cm")


# Plot top 10 for canonical structure and top 5 for other patterns
top_canonical_df <- subgenome_group_df[subgenome_group_df$Structure == "#----------------#####", ][1:10, ]
top_non_canonical_df <- plyr::ddply(
  subgenome_group_df[subgenome_group_df$Structure != "#----------------#####", ], 
  "Structure", function(block){
    return(block[1:5, ])
  })

top_df <- rbind(top_canonical_df, top_non_canonical_df)
top_df <- top_df[order(top_df$Structure, -top_df$JuncNum, decreasing = TRUE), ]
top_df <- na.omit(top_df)
plot_str_li <- build_str_df(top_df, chrom_size=chrom_size)
p <- ggplot() +
  geom_rect(data=gene_ref_df, aes(xmin=Start, xmax=End, ymin=-0.3, ymax=0.3), fill="black") +
  geom_text(data=gene_ref_df, aes(x=Mid, y=1, label=Name), color="black", size=1.2) +
  geom_rect(data=plot_str_li[["exon"]], aes(xmin=Start, xmax=End, ymin=-y-0.3, ymax=-y+0.3, fill=Str)) +
  geom_linerange(data=plot_str_li[["intron"]], aes(xmin=Start, xmax=End, y=-y, color=Str), size=0.3) +
  geom_rect(data=plot_str_li[["number"]], aes(xmin=num_start, xmax=num_start+1000*LogReadNum, ymin=-y-0.3, ymax=-y+0.3, fill=Str)) +
  geom_text(data=plot_str_li[["number"]], aes(x=num_start+1000*LogReadNum, y=-y, label=Label), color="black", size=1.2, hjust=-0.1) +
  geom_vline(xintercept = num_start, color="black", size=0.3) +
  scale_x_continuous(
    breaks = c(0, 15000, chrom_size, num_start+1000*log10(1+1), num_start+1000*log10(1+1000), num_start+1000*log10(1+1000000)), 
    labels = c("0", "15000", as.character(chrom_size), "1e0", "1e3", "1e6")
  ) +
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") +
  theme_bw() +
  theme(
    text = element_text(family="ArialMT", size=6),
    title = element_blank(),
    axis.text = element_text(color = "black"),
    axis.ticks.y = element_blank(),
    legend.key.size = unit(4, "mm"),
    legend.title = element_blank(),
    legend.text = element_text(family="ArialMT", size=6),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_line(color="black", size=0.3),
    axis.line.x = element_line(color="black", size=0.3),
    axis.text.y = element_blank(),
    axis.text.x = element_text(hjust=1)
  )
ggsave(file.path(oSG, "subgenome.str.top.pdf"), p, width = 12, height = 7.5, units = "cm")


