# Load ggplot2
library(ggplot2)
library(ggpubr)



barplot_fig <- function(fig, y_title){
  p1 <- ggbarplot(data, x = "Type", y = "Value",fill = 'Type',palette = c('#E7E6E6','#4472C4','#E7E6E6','#C00000'),
            add = c("mean_se"),
            width = 0.7)+ ylab(y_title)+xlab("")
  
  
  p2<-p1+geom_jitter(aes(Type,Value), shape = 21, size=3, color = "black",fill="black",width = 0.15)
   
  if(fig == 'S13')
    ggsave("C:/Users/jumaj/Documents/University/aneuploidy/Code for GITHUB/article_results/Figure S13/FigS13A.png", width = 7, height = 7)
  else if(fig == 'S15')
    ggsave("C:/Users/jumaj/Documents/University/aneuploidy/Code for GITHUB/article_results/Figure S15/FigS15A.png", width = 7, height = 7)
  else if(fig == 'S16')
  ggsave("C:/Users/jumaj/Documents/University/aneuploidy/Code for GITHUB/article_results/Figure S16/FigS16A.png", width = 7, height = 7)

}


data = read.csv('C:/Users/jumaj/Documents/University/aneuploidy/Code for GITHUB/article_results/Figure S13/FigS13A.csv')
barplot_fig('S13', 'Norm to GAPDH norm to siCTL')


data = read.csv('C:/Users/jumaj/Documents/University/aneuploidy/Code for GITHUB/article_results/Figure S15/FigS15A.csv')
barplot_fig('S15', 'Norm to GAPDH norm to siCTL')

data = read.csv('C:/Users/jumaj/Documents/University/aneuploidy/Code for GITHUB/article_results/Figure S16/FigS16A.csv')
barplot_fig('S16', 'Norm to GAPDH norm to siCTL')
