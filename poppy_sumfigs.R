#!/usr/bin/env Rscript

# Load command line arguments
args = commandArgs(trailingOnly=T)
infile = args[1]
outfile1 = args[2]
outfile2 = args[3]

# Load libraries
library(ggplot2)
library(RColorBrewer)

# Read input table
df = read.table(infile, header=T)

# Calculate colour intervals for MDF summary plot
if (sum(df$MDF > 0, na.rm=T) != 0){
  # If there are feasible pathways, calculate MDF intervals
  min_MDF = floor(min(df$MDF[df$MDF > 0], na.rm=T))
  max_MDF = floor(max(df$MDF[df$MDF > 0], na.rm=T))
  interval_seq = round(seq(min_MDF, max_MDF, length.out=5),1)
  df$interval = findInterval(df$MDF, interval_seq)
  # Create labels
  labels = c('Failed', '≤0',
             paste(
                   interval_seq[1:length(interval_seq)-1],
                   interval_seq[2:length(interval_seq)], sep='-'
             ),
             paste('>', interval_seq[length(interval_seq)], sep='')
  )
  interval = c(NA,0,1,2,3,4,5)
} else {
  # If there are no feasible pathways, specify intervals explicitly
  df$interval = ifelse(is.na(df$MDF), NA, 0)
  # Create labels
  labels = c('Failed', '≤0')
  interval = c(NA, 0)
}

# Create labels dataframe
label_df = data.frame(interval, factor(labels, rev(labels)))

# Add labels to dataframe
colnames(label_df) = c('interval', 'label')
df = merge(df, label_df)

# Prepare for plotting
df$x = rep('', nrow(df))

# Create the color scale
if (NA %in% df$MDF) {
  colors_vector = c('#525252','#E0E0E0')
} else {
  colors_vector = c('#E0E0E0')
}
# Add more colors if there are feasible pathways
if (sum(df$MDF > 0, na.rm=T) != 0){
  colors_vector = c(colors_vector, brewer.pal(5, 'OrRd'))
}

# Plot MDF summary
gp = ggplot(df, aes(x, fill=label, group=label))
gp = gp + geom_bar()
gp = gp + theme_bw()
gp = gp + scale_fill_manual(
  values=rev(colors_vector),
  guide=guide_legend(title="kJ/mol", reverse=F)
)
gp = gp + theme(
  axis.title.x=element_blank(),
  axis.text.x=element_blank(),
  axis.ticks.x=element_blank(),
  axis.title.y=element_blank()
)
ggsave(outfile1, gp, width=2.5, height=5)

# Plot length summary
df$length = factor(df$length, levels=sort.int(unique(df$length), decreasing=F))

gp = ggplot(df, aes(x, fill=length, group=length))
gp = gp + geom_bar()
gp = gp + theme_bw()
gp = gp + scale_fill_manual(
  values=rev(brewer.pal(max(length(unique(df$length)), 3), 'BuGn')),
  guide=guide_legend(title="Length\nin reactions", reverse=F)
)
gp = gp + theme(
  axis.title.x=element_blank(),
  axis.text.x=element_blank(),
  axis.ticks.x=element_blank(),
  axis.title.y=element_blank()
)
ggsave(outfile2, gp, width=2.5, height=5)
