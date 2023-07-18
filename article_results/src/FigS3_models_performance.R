library(data.table)
library(ggplot2)
library(reshape2)



get_performance_values <- function(data,performance_type){
  
  model_vector = vector()
  values_vector = vector()
  for(i in 1:10){
    model_vector <- c(model_vector, 'XGBoost')
    model_vector <- c(model_vector, 'Gradient boosting')
    model_vector <- c(model_vector, 'Random forest')
    model_vector <- c(model_vector, 'Bagging')
    model_vector <- c(model_vector, 'Logestic regression')
    
    if(performance_type == 'PRC'){
      values_vector <- c(values_vector, data$XGBoost.PR[i])
      values_vector <- c(values_vector, data$Gradient.Boosting.PR[i])
      values_vector <- c(values_vector, data$Random.Forest.PR[i])
      values_vector <- c(values_vector, data$Bagging.PR[i])
      values_vector <- c(values_vector, data$Logistic.Regression.PR[i])
    }
    else{
      values_vector <- c(values_vector, data$XGBoost.ROC[i])
      values_vector <- c(values_vector, data$Gradient.Boosting.ROC[i])
      values_vector <- c(values_vector, data$Random.Forest.ROC[i])
      values_vector <- c(values_vector, data$Bagging.ROC[i])
      values_vector <- c(values_vector, data$Logistic.Regression.ROC[i])
    }
  }
  
  return (list(model_vector,values_vector))
  
}


plot_boxplot <- function(data,performance_type,intercept_line){
  
  ggplot(data, aes(x=Type, y=values_vector, fill=Type)) +
    ylim(0,1)+
    xlab("")+
    ylab(performance_type)+
    ggtitle("")+
    theme( axis.title=element_text(size=20))+
    theme(axis.text.x=element_text(angle=90, hjust=1 , vjust=0.5 , size = 18, color = 'black'))+
    theme(axis.text.y=element_text(size = 18,color = 'black'))+
    theme(panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),panel.border = element_blank())+
    geom_boxplot(width=0.6, color='black',position = position_dodge(width=0.8))+
    scale_fill_manual(values=c('#66c2a5','#e78ac3','#fc8d62','#a6d854','#8da0cb'))+
    theme(legend.position = 'none',legend.text = element_text(size=18))+
    geom_hline(yintercept = intercept_line, size = 1, linetype = 2)
  
  save_file_path <- paste("../Figure S3/",model,sep='')
  if (performance_type == 'PRC'){
    save_file_path <- paste(save_file_path,'_PCR_boxplot.png',sep='')
  }
  else{
    save_file_path <- paste(save_file_path,'_ROC_boxplot.png',sep='')
  }
  ggsave(save_file_path)
}




models_list <- c('Loss versus Neutral','Loss versus Rest','Gain versus Neutral','Gain versus Rest')
expected_threshold <- c(0.42,0.33,0.32,0.21)
performance_type <- c('PRC','ROC')
index <- 1
for (model in models_list){
  model_file <- paste("../Figure S3/",model,sep = "")
  model_file <- paste(model_file,'__Ten_fold_results.csv',sep = "")
  data =read.csv(model_file)
  
  for (p_t in performance_type){
    per_values <- get_performance_values(data, p_t)
    model_vector <- unlist(per_values[1])
    values_vector <- unlist(per_values[2])
    df <- data.frame(model_vector, values_vector)
    df$Type<- factor(df$model_vector, c('XGBoost','Gradient boosting','Random forest','Bagging','Logestic regression'))
    
    if (p_t == 'PRC'){
      plot_boxplot(df,p_t,expected_threshold[index])
    }
    else{
      plot_boxplot(df,p_t,0.5)
    }
    
  }
  index <- index + 1
  
  
  
}
