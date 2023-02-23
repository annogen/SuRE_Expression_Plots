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
library(optparse)

# Parsing Opts
option_list = list(
  make_option(c("-f","--normfile"), type="character", 
              help="Path to the combined and subsetted normalised file [required]",default = FALSE),
  make_option(c("-o", "--outdir"), type="character", 
              help="Output directory [required]",default = FALSE),
  make_option(c("-p", "--position"), type="character", 
              help="Position of the mutation of interest [required]",default = FALSE),
  make_option(c("-c", "--chromosome"), type="character", 
              help="Chromosome of interest [required]",default = FALSE),
  make_option(c("-j", "--jitter"), type="character", default=1, 
              help="itter between the fragments"),
  make_option(c("-z", "--zerooffset"), type="character", default=1, 
              help="Ofset for fragments with zero expression, so that we can look at them seperately.")); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Input data
pos = as.numeric(opt$position)
chr = opt$chromosome
norm.file = opt$normfile


# Parameters  
jit = as.numeric(opt$jitter)
zero.exp.offset = as.numeric(opt$zerooffset)

# Output files
plotout = file.path(opt$outdir,paste0(chr,"_",pos,"_plot.png"))
statsout = file.path(opt$outdir,paste0(chr,"_",pos,"_stats.png"))


# Stats to be shown with the plot 
get_stats = function(chr,pos,snpvardf,statsout){ 
  
  # Split the dataframe by library and allele
  sp.snpvardf = split(snpvardf,list(snpvardf$allele,snpvardf$library))
  
  # Check if there is one library - this will affect the ouput format
  if(length(unique(snpvardf$library))>1){
    # We also want to display combined library information
    sp.snpvardf$`Alternate.Combined` = snpvardf[snpvardf$allele == "Alternate",]
    sp.snpvardf$`Reference.Combined` = snpvardf[snpvardf$allele == "Reference",]
    sp.snpvardf$`Reference.Combined`$library = "Combined"
    sp.snpvardf$`Alternate.Combined`$library = "Combined"
    
    # Calculate statistics 
    stats.sp.snpvardf = c(unlist(lapply(sp.snpvardf, function(x) as.character(dim(x)[1]))),
                          unlist(lapply(sp.snpvardf, function(y) as.character(round(mean(y$expression),digits = 2)))),
                          unlist(lapply(sp.snpvardf, function(y) as.character(round(100*(sum(as.numeric(y$expression) == 0)/dim(y)[1]))))))
    
    # Make a dataframe which is kable table compatible
    stats = data.frame("allele"=unlist(lapply(strsplit(names(stats.sp.snpvardf),"\\."), function(x) x[1])),
                       "Library"=unlist(lapply(strsplit(names(stats.sp.snpvardf),"\\."), function(x) x[2])),
                       "Value"=as.character(stats.sp.snpvardf),
                       "Statistic"=c(rep("Number of barcoded fragments overlapping the allele",(length(unique(snpvardf$library))+1)*2 ),
                                     rep("Mean expression of all barcoded fragments overlapping the allele",(length(unique(snpvardf$library))+1)*2 ),
                                     rep("Percentage of barcoded fragments with zero expression overlapping the allele",(length(unique(snpvardf$library))+1)*2 ))) 
  }else{
    # Calculate statistics 
    stats.sp.snpvardf = c(unlist(lapply(sp.snpvardf, function(x) as.character(dim(x)[1]))),
                          unlist(lapply(sp.snpvardf, function(y) as.character(round(mean(y$expression),digits = 2)))),
                          unlist(lapply(sp.snpvardf, function(y) as.character(round(100*(sum(as.numeric(y$expression) == 0)/dim(y)[1]))))))
    
    # Make a dataframe which is kable table compatible
    stats = data.frame("allele"=unlist(lapply(strsplit(names(stats.sp.snpvardf),"\\."), function(x) x[1])),
                       "Library"=unlist(lapply(strsplit(names(stats.sp.snpvardf),"\\."), function(x) x[2])),
                       "Value"=as.character(stats.sp.snpvardf),
                       "Statistic"=c(rep("Number of barcoded fragments overlapping the allele",length(unique(snpvardf$library))*2 ),
                                     rep("Mean expression of all barcoded fragments overlapping the allele",length(unique(snpvardf$library))*2 ),
                                     rep("Percentage of barcoded fragments with zero expression overlapping the allele",length(unique(snpvardf$library))*2 ))) 
    
  }
  
  # Define the levels , so that the order does not change when you use table or other fuctions 
  stats$Library = factor(stats$Library, levels =  unique(stats$Library))
  
  # Knit pretty table and save the table as PNG  
  s = kable(stats[c("allele","Value")],
            col.names = c("","")) %>% 
    kable_material(full_width = F) %>%
    pack_rows(index = table(fct_inorder(stats$Statistic))) %>%
    pack_rows(index = rep(table(split(stats,stats$Statistic)[[1]]$Library),length(split(stats,stats$Statistic)))) %>%
    kable_styling(font_size = 20)  %>%
    save_kable(file = statsout,zoom = 1.5)
  
  return()}

# Plot main plot
get_plot = function(chr,pos,snpvardf,offset,zero.exp.offset){
  # Position of mutation as you want it to be shown on the plot
  MUT = paste0("Position : ", chr,":",pos )
  
  # Min Y axix, all values have been shifted by offset and so is the min value but, 0 expression values have an additional chnage of zero.exp.offset in linear scale.
  # We set a Min Y for setting limits while plotting
  ymin = offset - zero.exp.offset
  # Maximum  ymax for setting limits while plotting
  ymax = ceiling(max(snpvardf$expression.jitter))
  
  # assign colours for each allele
  group.colors = c("Alternate" = "salmon", "Reference" = "darkturquoise")
  
  # Plot the fragment level expression plot
  g = ggplot(snpvardf, aes(x = start, xend = end , y  = expression.jitter , yend = expression.jitter, colour = as.factor(allele) ) ) + 
    geom_point(size = 0.01) + 
    geom_segment(aes(linetype = library )) + 
    geom_vline(xintercept = pos) + 
    xlab(MUT) +
    ylab("Normalised Expression") + 
    scale_color_manual(name="",values = group.colors ) + 
    scale_linetype_discrete(name="") +
    theme_classic() +
    # Y axis has been now transformed to log10 scale but the ticks are still in linear scale. This mean, the placement of the ticks is not linear.
    scale_y_continuous(trans='log10',
                       breaks = c(ymin, round(unique(10^(seq(log10(offset+1),log10(ymax), by = 0.4)))) ) , 
                       labels = c("0", round(unique(10^(seq(log10(offset+1),log10(ymax), by = 0.4))) - offset ) )) +
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
          legend.position = "top",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) + 
    ggtitle(label =  paste0("Normalised Expression plot for ",chr,":",pos),
            subtitle = "Fragment Level View")
  
  # Violin plot for refrence and alternate allel. This plot will help us see the difference between REF AND ALT when comparnig raQTLs to other mutations
  exp_vp = ggplot(snpvardf,aes(x = allele, y = (expression+1), colour = allele)) + geom_violin() +
    geom_sina() +
    theme_classic() + 
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) + 
    scale_colour_manual(name="",values = group.colors ) + 
    ylab("Normalised Expression") + xlab("")  +
    scale_y_continuous(trans='log10',
                       breaks = c(1,unique(round(10^(seq(log10(offset+1),log10(ymax), by = 0.4)))) - offset + 1) ,
                       labels = c(0,unique(round(10^(seq(log10(offset+1),log10(ymax), by = 0.4)))) - offset)  ) +
    stat_summary(fun = function(y) log10(mean(10 ^ y)),
                 geom = "point",
                 colour = "black") +
    ggtitle(label =  " ",
            subtitle = "Normalised expression")
  
  
  # Histogram scaled to 100% wrt to total counts in Reference and Alternate respectively.
  # This will help us see the biases in start sites of of alleles.
  start_vp = ggplot(snpvardf, aes(x=start , fill = allele)) + scale_colour_manual(name="",values = group.colors ) +
    geom_histogram(aes(y=100*(after_stat(count)/sum(after_stat(count)))) , colour = "black", bins = 50, alpha = 0.8) + 
    theme_classic() + 
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) + 
    ylab("Coverage (in %)") + xlab("") +
    ggtitle(label =  " ",
            subtitle = "Left most coordinate of fragments")
  
  final_plt = grid.arrange(arrangeGrob(g, nrow = 1, ncol = 1),
                           arrangeGrob(exp_vp,start_vp, ncol = 2 , nrow = 1))
  return()
}



# Combine the norm.subset.files
df = fread(norm.file)

# snpvar 0 = HAPLOTYPE 1 ( REF )
# snpvar 1 = HAPLOTYPE 2 ( AL )
# keep only haplotype 1 and 2, ie, keep only snpvar == 0 or 1
snpvardf = df[df$snpvar %in% c(0,1),]
# Add allele name, Reference or Alternate based on 0 and 1 , or Haplotype 1 and Haplotype 2
snpvardf$allele = ifelse(snpvardf$snpvar == 0, "Reference","Alternate")

# How many barcodes do you discard
system(paste0("echo About ",round(100*(sum(!unique(df$bc) %in% unique(snpvardf$bc))/length(unique(df$bc)))),"% of the barcoded fragments originally selected to be considered, did not have correct allele." ))

# zero.exp.offset is introduced to shift the non zero values further apart from the zero values ( expression values ) , to show the contrast.
# Jitter added so that we can see the fragments better
# Adding jitter. +-(factor * d/5)  where d is about the smallest difference between x values
snpvardf$expression.jitter = jitter(ifelse(snpvardf$expression==0,-zero.exp.offset,snpvardf$expression), factor = jit) 

# Making expression.jitter compatible for log scale
# Calculating an offset to make sure everything is positive and everything is compatible with log scale
# min of snpvardf$expression.jitter would be zero.exp.offset + jitter. This value would be negative
offset =  ceiling(abs(min(snpvardf$expression.jitter))) + zero.exp.offset

# By adding the abs of offset we shift the whole graph up by offset
snpvardf$expression.jitter = snpvardf$expression.jitter + abs(offset)

system(paste0("echo There is an offset added to expression of all fragments only for plotting purposes. The offset is ",round(abs(offset), digits = 3)))

# Calculate stats and save stats as PNG
get_stats(chr,pos,snpvardf,statsout)

# Plot the expression plot and save as PNG
# High Resolution PNG
png(filename = plotout, width = 200, height = 200, units = "mm", res = 300)
get_plot(chr,pos,snpvardf,offset,zero.exp.offset) 
dev.off()



