library(readr)
library(argparse)
library(dplyr)
library(reshape2)
library(argparse)
library(RColorBrewer)
library(ggplot2)
parser <- ArgumentParser()
parser$add_argument("--MS-data", nargs="*", required=TRUE, type="character", dest = "ms_data", metavar="ms_data.bed")
parser$add_argument("--sample", nargs="*", required=TRUE, type="character", dest = "sample", metavar="sample")
parser$add_argument("-o", "--output", required=TRUE, type="character", dest = "output", metavar="output")
args <- parser$parse_args(args)

load_bed <- function(fname, s){
  df <- read_delim(fname, "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE, col_types = cols(`X11` = col_character(), `X12` = col_character()))
  names(df) <- c("Chrom", "Start", "End", "Name", "Score", "Strand", "txStart", "txEnd", "color", "blockNum", "blockLen", "blockStart", "Tag")
  df <- df[,c("Start", "End", "Name", "Strand", "blockNum", "blockLen", "blockStart")]
  df$Sample <- s
  df <- df[df$blockNum==2,]
  df <- df[df$Start>100 & df$Start<20000 & df$End>20000, ]
  df <- df[df$Strand=="+",]
  df$Pipetide <- sapply(strsplit(df$Name, "[|]"), function(x){return(x[1])})
  df$PipetideNum <- sapply(strsplit(df$Name, "[|]"), function(x){return(as.integer(x[3]))})
  df$LeftBlockLen <- sapply(strsplit(df$blockLen, "[,]"), function(x){return(as.integer(x[1]))})
  df$RightBlockLen <- sapply(strsplit(df$blockLen, "[,]"), function(x){return(as.integer(x[2]))})
  df$LeftBlockEnd <- df$Start + df$LeftBlockLen
  df$RightBlockStart <- df$End - df$RightBlockLen
  return(df)
}
MS_bed_li <- list()
for(indx in 1:length(args$sample)){
MS_bed_li[[args$sample[indx]]] <- load_bed(args$ms_data[indx], args$sample[indx])
}
MS_bed_df <- do.call(rbind, MS_bed_li)
MS_bed_df <- MS_bed_df[MS_bed_df$LeftBlockLen>10 & MS_bed_df$RightBlockLen>10,]
MS_js_df <- MS_bed_df %>% group_by(Sample, LeftBlockEnd, RightBlockStart) %>% summarise(PipeNum=sum(PipetideNum))

p <- ggplot(MS_js_df, aes(x=LeftBlockEnd, y=0, xend=RightBlockStart, yend=0)) +
  geom_curve(color="black", curvature=-0.2, size=0.1) +
  facet_grid(Sample~.) +
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

MS_js_data <- MS_bed_df %>% group_by(Pipetide, Start, End, blockLen, blockStart) %>% summarise(SampleNum=n())
enriched_MS_js_data <- MS_js_data[MS_js_data$SampleNum>1,]
enriched_MS_bed_df <- inner_join(enriched_MS_js_data, MS_bed_df)
## Check peptides with both b ion and y ion evidence support
write_tsv(enriched_MS_bed_df, file.path(args$output, "ORF1ab.MS.tsv"))

