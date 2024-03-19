library(ggplot2)
library(corrplot)


### Fig. S3 panel A - heatmap correlation
setwd("C:/Users/jumaj/Documents/University/aneuploidy/Code for GITHUB/src/Supplementary analyses/Fig S3/")


data=read.csv("feature_correlation_results.csv", header = T, sep = ',')
rownames(data)<-unlist(data$feature)
data$feature<-NULL

corrplot(as.matrix(data), method = 'ellipse', order = 'AOE', type = 'upper',tl.cex = 1, tl.col = "black", col = COL2('BrBG', 200))

corrplot(as.matrix(data), method = 'number', order = 'AOE', type = 'upper',tl.cex = 1, number.cex = 0.8,tl.col = "black", col = COL2('BrBG', 200))




### Fig. S3 panel C - features PCA




library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)
library(factoextra)
library(ggfortify)
library(plotly)

setwd("C:/Users/jumaj/Documents/University/aneuploidy/Code for GITHUB/article_results/Data/")
data=read.csv("Model_dataset.csv", header = T, sep = ',')

data$Arm<-NULL
data$Type<-NULL
data$Label<-NULL
data <- data %>% mutate_if(is.numeric, function(data) ifelse(is.na(data), median(data, na.rm = T), data))

# calculate PCA
res.pca<-prcomp(data, scale = FALSE)

# get first three components
components <- res.pca[["x"]]
components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3
components$PC1
data=read.csv("Model_dataset.csv", header = T, sep = ',')



components = cbind(components, data$Label)


fig <- plot_ly(components, x = ~PC1, y = ~PC2, z = ~PC3, color = ~as.character(data$Label), colors = c('blue','#F5F5DC','red'), mode = 'markers')%>%
  add_trace(marker = list(line = list(color = 'black',width = 1)))

fig <- fig %>%
  layout(
    scene = list(bgcolor = "white")
  )

fig

