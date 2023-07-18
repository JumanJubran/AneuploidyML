library(ggplot2)
library(ggsignif)
library(plyr)


bar_plotting <- function(data,data_type){
  
  data <- ddply(data, .(Gained),transform, pos = (cumsum(Value) - (0.5 * Value))/100)
  ggplot(data, aes(factor(Gained, level = c('No Gain','Gain')), y = Value, fill = Lost)) + 
    
    geom_bar(position = "fill",stat = "identity", width = 0.7, colour = 'black', size = 0.8) +
    ylab("Percentage of 4p loss")+
    xlab("")+
    scale_fill_manual(values=c('#4472C4','#E3CFCF'))+
    theme_classic()+
    theme(axis.title=element_text(size=20))+
    theme(axis.text.x=element_text(angle=0, hjust=0.5 , vjust=0 , size = 18, color = 'black'))+
    theme(axis.text.y=element_text(size = 18,color = 'black'))+
    theme(legend.position = 'none')+
    geom_signif(y_position = c(1.05), xmin = c(1), xmax = c(2),annotation = c("****"), tip_length = 0.0001, size = 1, textsize = 7)+
    geom_text(aes(label = Per), size = 5, y = data$pos, position = "identity")
  
  
  if (data_type == 'TCGA'){
    ggsave('../Figure 5/Fig5G_TCGA_cooccurence_barplot.png')
  }
  else{
    ggsave('../Figure 5/Fig5H_DepMap_cooccurence_barplot.png')
  }
  

}


data=read.csv("Fig5G_TCGA_cooccurence.csv", header = T, sep = ',')
bar_plotting(data,'TCGA')
data=read.csv("Fig5H_DepMap_cooccurence.csv", header = T, sep = ',')
bar_plotting(data,'DepMap')
