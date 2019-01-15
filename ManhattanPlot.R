# Plot SNPS
# Libraries ====
library(readr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(magrittr)
library(tidyr)

# Bed file should be in the following order (no header)
# chr start end region_id counts

file <- '~/Desktop/plots/36087_met_somatic_Mutect2.windows.bed'
x_label = 'Chromosomes'
y_label = 'Count Density'

output <- sub('.bed','.png',basename(file))
title <- sub('.png','',output)

snps<-read.table(file,sep="\t",header=F)
colnames(snps)<-c('chr','binstart','binend','binid','count')
# remove chr prefix from chr column first
snps$chr <- sub('^chr','',snps$chr)
# put the chromosomes in the good order: chr1, chr2, chr22, chrX
goodChrOrder <- paste("",c(1:22,"X","Y"),sep="")
snps$chr <- factor(snps$chr,levels=goodChrOrder)
snps <- snps[!is.na(names(snps))]

snps %<>% drop_na()

gg.manhattan <- function(df, col, title, x_lab = x_label, y_lab = y_label){
  df.tmp <- df %>% 
    group_by(chr) %>% 
    summarise(chr_len=max(binend)) %>% 
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    left_join(df, ., by=c("chr"="chr")) %>%
    arrange(chr, binstart) %>%
    mutate( BPcum=binstart+tot) %>%
    mutate(log2count = log2(count))
  
  axisdf <- df.tmp %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2, separation.line=min(BPcum))
  shadingsdf <- df.tmp %>% group_by(chr) %>% summarise(min = min(BPcum), max=max(BPcum))
  shadingsdf$isShaded <- c(0,1)
  shadingsdf <- filter(shadingsdf, isShaded == 1)
  
  ggplot(df.tmp, aes(x=BPcum, y=count)) +
    # custom X axis:
    scale_x_continuous(expand=c(0,0), label = axisdf$chr, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0,NA)) + # expand=c(0,0)removes space between plot area and x axis 
    
    # add plot and axis titles
    ggtitle(paste0(title)) +
    labs(x = x_lab,
         y = y_lab) +
    
    # separate chromosomes with shading
    geom_rect(data = shadingsdf,
              aes(xmin = min, xmax = max, ymin = -Inf, ymax = Inf),
              alpha = 0.4,
              fill='grey',
              inherit.aes = FALSE)+
    
    # Show all points
    geom_point(aes(color=as.factor(chr)), alpha=0.8, size=0.5) +
    scale_color_manual(values = rep(col, 22 )) +
    
  # Custom the theme:
  theme_bw(base_size = 10) + 
    theme( 
      plot.title = element_text(hjust = 0.5),
      legend.position="none",
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black"))
}

mypalette <- c("#E2709A", "#CB4577", "#BD215B", "#970F42", "#75002B") # chr color palette
ggsave(
  output,
  gg.manhattan(snps, col=mypalette, title=title),
  width = 10,
  height = 4,
  dpi = 1200
)


