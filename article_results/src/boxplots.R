library(ggplot2)
library(ggsignif)
library(plyr)


box_plotting <- function(asterisk_position, y_lim_low, y_lim_high,step_by, y_title, figure, significance,groups,colors_c){
  ggplot(data, aes(factor(type, level = groups),  value , fill=type))+ 
    geom_boxplot(width = 0.5, size = 0.6)+ 
    geom_point(position = position_dodge2(width = 0.2),alpha =0.5) +
    scale_y_continuous(breaks = seq(y_lim_low, y_lim_high, by = step_by), limits = c(y_lim_low, y_lim_high))+
    scale_fill_manual(values=colors_c)+
    
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
  else if(figure=='4J'){
    ggsave('C:/Users/jumaj/Documents/University/aneuploidy/Code for GITHUB/article_results/Figure 4/Fig4J_siKLF5_DLD1.png')
  }
  else if(figure=='S12_wt'){
    ggsave('C:/Users/jumaj/Documents/University/aneuploidy/Code for GITHUB/article_results/Figure S12/S12_wt.png')
  }
  else if(figure=='S12_Ts13'){
    ggsave('C:/Users/jumaj/Documents/University/aneuploidy/Code for GITHUB/article_results/Figure S12/S12_Ts13.png')
  }
  else if(figure=='S14C'){
    ggsave('C:/Users/jumaj/Documents/University/aneuploidy/Code for GITHUB/article_results/Figure S14/S14C.png')
  }
  else if(figure=='S15D'){
    ggsave('C:/Users/jumaj/Documents/University/aneuploidy/Code for GITHUB/article_results/Figure S15/S15D.png')
  }
  else if(figure=='S16D'){
    ggsave('C:/Users/jumaj/Documents/University/aneuploidy/Code for GITHUB/article_results/Figure S16/S16D.png')
  }
  
  
  
}

### Figure 4
groups<-c('DLD1 WT','DLD1 Ts13')
colors_c <- c('#C00000','#E7E6E6')

data = read.csv('C:/Users/jumaj/Documents/University/aneuploidy/Code for GITHUB/article_results/Figure 4/Fig4H_KLF5_mRNA.csv')
box_plotting(2.6, 0, 3,1, 'KLF5 mRNA levels', '4H', '**',groups,colors_c)

data = read.csv('C:/Users/jumaj/Documents/University/aneuploidy/Code for GITHUB/article_results/Figure 4/Fig4J_siKLF5_DLD1.csv')
box_plotting(1, 0, 1,0.25, 'Relative viability\nFollowing siKLF5', '4J', '*',groups,colors_c)


### Figure S12
groups<-c('DLD1 WT siCTL','DLD1 WT siKLF5')
colors_c <- c('#E7E6E6','#4472C4')
data = read.csv('C:/Users/jumaj/Documents/University/aneuploidy/Code for GITHUB/article_results/Figure S12/FigS12_WT.csv')
box_plotting(1.3, 0, 1.3,0.25, 'KLF5 mRNA levels', 'S12_wt', '**',groups,colors_c)

groups<-c('DLD1 Ts13 siCTL','DLD1 Ts13 siKLF5')
colors_c <- c('#E7E6E6','#C00000')
data = read.csv('C:/Users/jumaj/Documents/University/aneuploidy/Code for GITHUB/article_results/Figure S12/FigS12_Ts13.csv')
box_plotting(1.3, 0, 1.3,0.25, 'KLF5 mRNA levels', 'S12_Ts13', '*',groups,colors_c)


### Figure S14
groups<-c('WT siNEK3','DLD1 Ts13 siCTL')
colors_c <- c('#C00000','#4472C4')
data = read.csv('C:/Users/jumaj/Documents/University/aneuploidy/Code for GITHUB/article_results/Figure S14/FigS14C.csv')
box_plotting(1.3, 0, 1.3,0.25, 'Relative viability\nFollowing siCTL', 'S14C', 'n.s.',groups,colors_c)

### Figure S15
groups<-c('WT siNEK3','Ts13 siNEK3')
colors_c <- c('#C00000','#4472C4')
data = read.csv('C:/Users/jumaj/Documents/University/aneuploidy/Code for GITHUB/article_results/Figure S15/FigS15D.csv')
box_plotting(1.3, 0, 1.5,0.5, 'Relative viability\nFollowing siNEK3', 'S15D', 'n.s.',groups,colors_c)

### Figure S16
groups<-c('WT siTTC7A','Ts13 siTTC7A')
colors_c <- c('#C00000','#4472C4')
data = read.csv('C:/Users/jumaj/Documents/University/aneuploidy/Code for GITHUB/article_results/Figure S16/FigS16D.csv')
box_plotting(1.3, 0, 1.5,0.25, 'Relative viability\nFollowing siTTC7A', 'S16D', 'n.s.',groups,colors_c)




