# AUTHOR / DATE
# Vartika Bisht; January 31, 2023

## Plot expression plots and statistics using the subsetted SuRE counts data.
## These plots are saved as PNGs


# required libraries
library(data.table)
library(ggplot2)
library(kableExtra)
library(forcats)
library(gridExtra)
library(reshape2)
library(ggExtra)
library(ggforce)

# Input data
cols = snakemake@params[["cols"]]
libs = snakemake@params[["libs"]]
mut = snakemake@params[["mut"]]
libname = snakemake@params[["libname"]]
jit = as.numeric(snakemake@params[["jitter"]])
zero.exp.offset = as.numeric(snakemake@params[["zero_exp_offset"]])
norm.subset.files = snakemake@input
plotout = snakemake@output[["plotout"]]
statsout = snakemake@output[["statsout"]]


# Reformat data frame for plotting
reform = function(f,cols,libname){
  df = fread(f, header = TRUE)
  result_df = df[,.SD, .SDcols = as.character(cols)]
  colnames(result_df) = names(cols)
  if(dim(df)[1] != 0){
    result_df$ipcr.norm.sum = rowSums(df[,.SD, .SDcols = grep("ipcr.norm.sum",colnames(df))])/2
    result_df$Allele = df[,.SD,.SDcols = as.character(cols["MUTVAR"])]
    result_df = result_df[result_df$Allele %in% c(0,1),]
    result_df$Allele = ifelse(result_df$Allele == 0, "Reference Allele" , "Alternate Allele")
    result_df$lib = libname[unique(result_df$lib)] }
  return(result_df)}

# Stats to be shown with the plot 
get_stats = function(mut,file.list,statsout){ 
  MUT = paste0("Position : ", strsplit(mut,"_")[[1]][1],":",as.numeric(strsplit(mut,"_")[[1]][2]))
  sp.file.list = split(file.list,list(file.list$Allele,file.list$lib))

  # Check if there is one library - this will affect the ouput format
  if(length(unique(file.list$lib))>1){
  sp.file.list$`Alternate Allele.Combined` = file.list[file.list$Allele == "Alternate Allele",]
  sp.file.list$`Reference Allele.Combined` = file.list[file.list$Allele == "Reference Allele",]
  sp.file.list$`Reference Allele.Combined`$lib == "Combined"
  sp.file.list$`Alternate Allele.Combined`$lib == "Combined"
  
  stats.sp.file.list = c(unlist(lapply(sp.file.list, function(x) as.character(dim(x)[1]))),
                             unlist(lapply(sp.file.list, function(y) as.character(round(mean(y$ipcr.norm.sum),digits = 2)))),
                             unlist(lapply(sp.file.list, function(y) as.character(round(100*(sum(as.numeric(y$ipcr.norm.sum) == 0)/dim(y)[1]))))))
  
  stats = data.frame("Allele"=unlist(lapply(strsplit(names(stats.sp.file.list),"\\."), function(x) x[1])),
                     "Library"=unlist(lapply(strsplit(names(stats.sp.file.list),"\\."), function(x) x[2])),
                     "Value"=as.character(stats.sp.file.list),
                     "Statistic"=c(rep("Number of barcoded fragments overlapping the allele",(length(unique(file.list$lib))+1)*2 ),
                                   rep("Mean expression of all barcoded fragments overlapping the allele",(length(unique(file.list$lib))+1)*2 ),
                                       rep("Percentage of barcoded fragments with zero expression overlapping the allele",(length(unique(file.list$lib))+1)*2 ))) 
  }else{
    stats.sp.file.list = c(unlist(lapply(sp.file.list, function(x) as.character(dim(x)[1]))),
                           unlist(lapply(sp.file.list, function(y) as.character(round(mean(y$ipcr.norm.sum),digits = 2)))),
                           unlist(lapply(sp.file.list, function(y) as.character(round(100*(sum(as.numeric(y$ipcr.norm.sum) == 0)/dim(y)[1]))))))
    
    stats = data.frame("Allele"=unlist(lapply(strsplit(names(stats.sp.file.list),"\\."), function(x) x[1])),
                       "Library"=unlist(lapply(strsplit(names(stats.sp.file.list),"\\."), function(x) x[2])),
                       "Value"=as.character(stats.sp.file.list),
                       "Statistic"=c(rep("Number of barcoded fragments overlapping the allele",length(unique(file.list$lib))*2 ),
                                     rep("Mean expression of all barcoded fragments overlapping the allele",length(unique(file.list$lib))*2 ),
                                     rep("Percentage of barcoded fragments with zero expression overlapping the allele",length(unique(file.list$lib))*2 ))) 
    
  }

  
  # Knit pretty table and save the table as PNG  
  s = kable(stats[c("Allele","Value")],
            caption = MUT,
            col.names = c("","")) %>% 
    kable_material(full_width = F) %>%
    pack_rows(index = table(fct_inorder(stats$Statistic))) %>%
    pack_rows(index = rep(table(split(stats,stats$Statistic)[[1]]$Library),length(split(stats,stats$Statistic)))) %>%
    kable_styling(font_size = 20)  %>%
    save_kable(file = statsout,zoom = 1.5)

  return()}

# Plot main plot
get_plot = function(mut,file.list,offset,zero.exp.offset){
  # Position of mutation as you want it to be shown on the plot
  MUT = paste0("Position : ", strsplit(mut,"_")[[1]][1],":",as.numeric(strsplit(mut,"_")[[1]][2]))

  # Min Y axix, all values have been shifted by offset and so is the min value but, 0 expression values have an additional chnage of zero.exp.offset in linear scale.
  # We set a Min Y for setting limits while plotting
  ymin = offset - zero.exp.offset
  # Maximum  ymax for setting limits while plotting
  ymax = ceiling(max(file.list$ipcr.norm.sum.jitter))
  
  # assign colours for each allele
  group.colors = c("Alternate Allele" = "salmon", "Reference Allele" = "darkturquoise")
  
  # Plot the fragment level expression plot
  g = ggplot(file.list, aes(x = start, xend = end , y  = ipcr.norm.sum.jitter , yend = ipcr.norm.sum.jitter, colour = as.factor(Allele) ) ) + 
    geom_point(size = 0.01) + 
    geom_segment(aes(linetype = lib )) + 
    geom_vline(xintercept = as.numeric(strsplit(mut,"_")[[1]][2])) + 
    xlab(MUT) +
    ylab("Normalised Expression") + 
    scale_color_manual(name="",values = group.colors ) + 
    scale_linetype_discrete(name="") +
    theme_classic() +
    # Y axis has been now transformed to log10 scale but the ticks are still in linear scale. This mean, the placement of the ticks is not linear.
    scale_y_continuous(trans='log10',
                       breaks = c(ymin, unique(round(10^(seq(log10(offset+1),log10(ymax), by = 0.2)))) ) , 
                       labels = c("0", unique(round(10^(seq(log10(offset+1),log10(ymax), by = 0.2)))) - offset ) ) +
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
          legend.position = "top",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) + 
    ggtitle(label =  paste0("Normalised Expression plot for Chr",strsplit(mut,"_")[[1]][1],":",as.numeric(strsplit(mut,"_")[[1]][2])),
            subtitle = "Fragment Level View")
  
  # Violin plot for refrence and alternate allel. This plot will help us see the difference between REF AND ALT when comparnig raQTLs to other mutations
  exp_vp = ggplot(file.list,aes(x = Allele, y = (ipcr.norm.sum+1), colour = Allele)) + geom_violin() +
    geom_sina() +
    theme_classic() + 
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) + 
    scale_colour_manual(name="",values = group.colors ) + 
    ylab("Normalised Expression") + xlab("")  +
    scale_y_continuous(trans='log10',
                       breaks = c(1,unique(round(10^(seq(log10(offset+1),log10(ymax), by = 0.2)))) - offset + 1) ,
                       labels = c(0,unique(round(10^(seq(log10(offset+1),log10(ymax), by = 0.2)))) - offset)  ) +
    stat_summary(fun = function(y) log10(mean(10 ^ y)),
                 geom = "point",
                 colour = "black") +
    ggtitle(label =  "Combined View",
            subtitle = "Normalised expression between alleles")
  
  
  # Violin plot for library level 
  exp_vp_lib = ggplot(file.list,aes(x = Allele, y = (ipcr.norm.sum+1), colour = Allele)) + geom_violin() +
    geom_sina() +
    theme_classic() + 
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) + 
    scale_colour_manual(name="",values = group.colors ) + 
    ylab("Normalised Expression") + xlab("")  +
    scale_y_continuous(trans='log10',
                       breaks = c(1,unique(round(10^(seq(log10(offset+1),log10(ymax), by = 0.2)))) - offset + 1) ,
                       labels = c(0,unique(round(10^(seq(log10(offset+1),log10(ymax), by = 0.2)))) - offset)  ) +
    stat_summary(fun = function(y) log10(mean(10 ^ y)),
                 geom = "point",
                 colour = "black") +
    ggtitle(label =  "Combined View",
            subtitle = "Normalised expression between alleles") + facet_wrap(~lib)
  
  # Histogram scaled to 100% wrt to total counts in reference allele and alternate allele respectively.
  # This will help us see the biases in start sites of of alleles.
  start_vp = ggplot(file.list, aes(x=start , fill = Allele)) + scale_colour_manual(name="",values = group.colors ) +
    geom_histogram(aes(y=100*(after_stat(count)/sum(after_stat(count)))) , colour = "black", bins = 50, alpha = 0.8) + 
    theme_classic() + 
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) + 
    ylab("Propotion of start sites (in %)") + xlab("") +
    ggtitle(label =  "Combined View",
            subtitle = "Start-sites of fragment between alleles")
  
  start_vp_lib = ggplot(file.list, aes(x=start , fill = Allele)) + scale_colour_manual(name="",values = group.colors ) + 
    facet_wrap(~lib) +
    geom_histogram(aes(y=100*after_stat(count/tapply(count, PANEL, sum)[PANEL])) , colour = "black", bins = 50, alpha = 0.8) + 
    theme_classic() + 
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) + 
    ylab("Propotion of start sites (in %)") + xlab("") +
    ggtitle(label =  "Library Level View",
            subtitle = "Start-sites of fragment between alleles")
  
  return(grid.arrange(arrangeGrob(g, nrow = 1, ncol = 1),
                      arrangeGrob(exp_vp,exp_vp_lib, ncol = 2 , nrow = 1),
                      arrangeGrob(start_vp,start_vp_lib, ncol = 2 , nrow = 1)))
}



# Combine the norm.subset.files
file.list = list()
for(f in norm.subset.files){
  file.list = c(file.list,list(reform(f,cols,libname)))}
file.list = Filter(function(x) dim(x)[1] > 0, file.list)
file.list = rbindlist(file.list)

# zero.exp.offset is introduced to shift the non zero values further apart from the zero values ( expression values ) , to show the contrast.
# Jitter added so that we can see the fragments better
# Adding jitter. +-(factor * d/5)  where d is about the smallest difference between x values
file.list$ipcr.norm.sum.jitter = jitter(ifelse(file.list$ipcr.norm.sum==0,-zero.exp.offset,file.list$ipcr.norm.sum), factor = jit) 
# Making ipcr.norm.sum.jitter compatible for log scale
# Calculating an offset to make sure everything is positive and everything is compatible with log scale
offset =  abs(round(min(file.list$ipcr.norm.sum.jitter))) + 1
file.list$ipcr.norm.sum.jitter = file.list$ipcr.norm.sum.jitter + offset

# Calculate stats and save stats as PNG
get_stats(mut,file.list,statsout)
# Plot the expression plot and save as PNG
# High Resolution PNG
png(filename = plotout, width = 350, height = 400, units='mm', res = 300)
get_plot(mut,file.list,offset,zero.exp.offset) 
dev.off()
