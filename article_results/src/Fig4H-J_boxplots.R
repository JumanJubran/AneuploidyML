library(ggplot2)
library(ggsignif)
library(plyr)


box_plotting <- function(asterisk_position, y_lim_low, y_lim_high,step_by, y_title, figure, significance){
  ggplot(data, aes(factor(type, level = c('DLD1 WT','DLD1 Ts13')),  value , fill=type))+ 
    geom_boxplot(width = 0.5, size = 0.6)+ 
    geom_point(position = position_dodge2(width = 0.2),alpha =0.5) +
    scale_y_continuous(breaks = seq(y_lim_low, y_lim_high, by = step_by))+
    scale_fill_manual(values=c('#C00000','#D8D6D6'))+
    
    # Remove background color
    theme_classic()+
    
    # Remove axis labels and ticks, in addition to legends
    ylab(y_title)+
    xlab("")+
    theme( axis.title=element_text(size=20))+
    theme(axis.text.x=element_text(angle=0, hjust=0.5 , vjust=0 , size = 18, color = 'black'))+
    theme(axis.text.y=element_text(size = 18,color = 'black'))+
    theme(legend.position = 'none')+
    
    # Add the brackets of the significance
    geom_signif(
      y_position = c(asterisk_position), xmin = c(1), xmax = c(2),
      annotation = c(significance), tip_length = 0.05, size = 1, textsize = 8
    )
  if(figure=='4H'){
    ggsave('C:/Users/jumaj/Documents/University/aneuploidy/Code for GITHUB/article_results/Figure 4/Fig4H_KLF5_mRNA_levels.png')
  }
  else{
    ggsave('C:/Users/jumaj/Documents/University/aneuploidy/Code for GITHUB/article_results/Figure 4/Fig4J_siKLF5_DLD1.png')
  }
  
}


data = read.csv('../Figure 4/Fig4H_KLF5_mRNA.csv')
box_plotting(2.6, 1, 3,0.5, 'KLF5 mRNA levels', '4H', '**')

data = read.csv('../Figure 4/Fig4J_siKLF5_DLD1.csv')
box_plotting(1, 0.6, 1,0.1, 'Relative viability\nFollowing siKLF5', '4J', '*')

