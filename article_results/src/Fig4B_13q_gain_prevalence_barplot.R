library(ggplot2)
library(ggsignif)
library(plyr)

data=read.csv("C:/Users/jumaj/Documents/University/aneuploidy/Code for GITHUB/article_results/Figure 4/Fig4B_13q_aneuploidy_prevalence.csv", header = T, sep = ',')


data <- ddply(data, .(cell), transform, percent = value/sum(value) * 100)
data <- ddply(data, .(cell),transform, cumsum_ = cumsum(percent))
data <- ddply(data, .(cell),transform, pos = (cumsum_ - percent/2)/100)


ggplot(data, aes(x = factor(cell, levels = c('Other','CRC')), y = value, fill = factor(type, levels = c('Gain','Neutral','Loss')))) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.25))+
  geom_bar(position = "fill",stat = "identity", width = 0.7, colour = 'black', size = 0.8)+
  ylab("Relative arm score (%)")+
  xlab("")+
  scale_fill_manual(values=c('#C00000','#E7E6E6','#4472C4'))+
  theme_classic()+
  
  
  theme( axis.title=element_text(size=20))+
  theme(axis.text.x=element_text(angle=0, hjust=0.5 , vjust=0 , size = 18, color = 'black'))+
  theme(axis.text.y=element_text(size = 18,color = 'black'))+
  theme(legend.position = 'none')+
  geom_signif(y_position = c(1.1), xmin = c(1), xmax = c(2),annotation = c("****"), tip_length = 0.00005, size = 1, textsize = 7)+
  geom_text(aes(label = Per), size = 5, y = data$pos, position = "identity")

ggsave('C:/Users/jumaj/Documents/University/aneuploidy/Code for GITHUB/article_results/Figure 4/Fig4B_13q_aneuploidy_prevalence.png')
