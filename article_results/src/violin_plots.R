library(ggplot2)
library(ggsignif)
library(plyr)
library(ggbeeswarm)

violin_plot <- function(path, data, y_label, significance, figure, asterisk_position, jitt_fun,color, group1, group2,adjust_parameter, weidth_parameter){

  colnames(data)[2] <- "value"
  colnames(data)[6] <- "type"
  
  data["type"][data["type"] == "out group"] <- group1
  
  data["type"][data["type"] == "in group"] <- group2
  
  # PLOT
  ggplot(data, aes(factor(type, level = c(group1,group2)), value , fill=type)) + 
    
    geom_violin(trim = FALSE, adjust = adjust_parameter, width = weidth_parameter,colour = c('#E7E6E6'))+
    jitt_fun+
    
    scale_fill_manual(values=color)+
    
    # Remove background color
    theme_classic()+
    
    # Remove axis labels and ticks, in addition to legends
    ylab(y_label)+
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
  
    #Save file
    if (figure == '5D'){
      ggsave('../Figure 5/Fig5D_UCHL1_essentiality.png')
    }
    else if (figure == '5E'){
      ggsave('../Figure 5/Fig5E_UCHL3_expression.png')
    }
    else if (figure == '4D'){
      ggsave('../Figure 4/Fig4D_KLF5_essentiality.png')
    }
    else if (figure == '4F'){
      ggsave('../Figure 4/Fig4F_KLF5_expression.png')
    }
}


  
jitt_fun <-  geom_point(position = position_jitter(width = 0.04),alpha =0.5)  
color <- c('#C00000','#D8D6D6')

data=read.csv("../Figure 5/Fig5D_UCHL1 Gene Effect (Chronos) CRISPR (DepMap Public 22Q4+Score Chronos).csv", header = T, sep = ',')
violin_plot(path, data,"UCHL1 essentiality\n (CRISPR chrons score)", '***', '5D',0.6,jitt_fun,color,'CRC chr13q WT','CRC chr13q gain',0.5,0.7)

data=read.csv("../Figure 5/Fig5E_UCHL3 log2(TPM+1) Expression Public 22Q4.csv", header = T, sep = ',')
violin_plot(path, data,'Expression\n UCHL3 log2(TPM+1)', '****', '5E',10,jitt_fun,color,'CRC chr13q WT','CRC chr13q gain',0.5,0.7)



jitt_fun <- geom_beeswarm(alpha = 0.2, cex = 0.5, shape = 21)
color <- c('#007A78','#FFC745')

data=read.csv("../Figure 4/Fig4D_KLF5 Gene Effect (Chronos) CRISPR (DepMap Public 22Q4+Score Chronos).csv", header = T, sep = ',')
violin_plot(path, data,'KLF5 essentiality\n (CRISPR chrons score)', '****', '4D',1,jitt_fun,color,'Other','CRC',0.7, 1)

data=read.csv("../Figure 4/Fig4F_KLF5 log2(TPM+1) Expression Public 22Q4.csv", header = T, sep = ',')
violin_plot(path, data,'Expression\n log2(TPM+1)', '****', '4F',10,jitt_fun,color,'Other','CRC',0.7,1)

