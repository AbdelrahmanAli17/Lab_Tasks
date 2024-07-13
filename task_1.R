library(openxlsx)
library(dplyr)
library("pheatmap")
library("RColorBrewer")
library(ggplot2)

plasma <- read.xlsx("QMDiab_metabolomics_Preprocessed.xlsx" , sheet = "plasma")

plsma_annotation <- read.xlsx("QMDiab_metabolomics_Preprocessed.xlsx" , sheet = "plasma annotations")


#---------------------------Heatmap for the plasma data--------------------------------- 

data_heatmap <- plasma
data_heatmap <- data_heatmap[order(data_heatmap$T2D),]

#write.csv(data_heatmap , "data_fake.csv")

rownames(data_heatmap)<-plasma$`QMDiab-ID`

T2D <- data.frame( T2D = ifelse(grepl("1",data_heatmap$T2D),
                                    "with","without"),
                     row.names = row.names(data_heatmap))

data_heatmap <- data_heatmap[,-c(1 ,c(501:507))]

#setwd("D:/CV/Krumsiek_task/plots/")

pheatmap(as.matrix(data_heatmap)  , 
         scale = "column", 
         clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation" , 
         clustering_method = "ward.D", #to minimize the variance between clusters 
         fontsize = 15,  color =  rev(colorRampPalette(brewer.pal(11, "RdGy"))(256)),
         angle_col = 45 , na_col = "grey",
         cluster_cols = T ,
         cluster_rows = T,
         filename = "myheatmap_2.png",
         annotation_row = T2D,
         annotation_names_row  = F,
         show_rownames = F , show_colnames = F
)

#shapiro test 
means <- data.frame(apply(data_heatmap, 1, mean))
shapiro.test(means$apply.data_heatmap..1..mean.)


data(means)
png("histo.png")
hist(means$apply.data_heatmap..1..mean. ,main = "Histogram of plsma data" ) 
dev.off()
#why 

data(means)
png("qqplot.png")
plot(density(means$apply.data_heatmap..1..mean.) , )
dev.off()

#why 
library(ggpubr)
ggsave( "qq.png",
        ggqqplot(means$apply.data_heatmap..1..mean.) , height = 10 , width = 10)
#linear regression 

plasma_reg <- plasma
rownames(plasma_reg)<-plasma_reg$`QMDiab-ID`
plasma_reg<- plasma_reg[,-c(1)]

p_values <-rep(NA, ncol(plasma_reg[,-c(501:507)]))

for(i in 1:ncol(plasma_reg[,-c(501:507)])){
  met <- plasma_reg[,i]

  model <- lm(met ~ T2D + AGE + GENDER + BMI 
            , data =plasma_reg ) 
  p_values[i] <- summary(model)$coefficients[2, 4]

}

#plot(p_values)

significants <-cbind(data.frame(t(plasma_reg[,-c(501:507)])) , p_values) 

significants_pvalue <- significants %>% dplyr::select(p_values)


significants_pvalue$color <- ifelse(significants_pvalue$p_values < 0.05 , "Significant", "NON")


x <- ggplot(significants_pvalue , aes(y = p_values , x = rownames(significants_pvalue), colour = color)) + 
  
  geom_point(size = 3, alpha = 0.7) + 
  
  theme(legend.background = element_rect( linetype="solid", 
                                         colour ="black")) + 
  theme( axis.line.y = element_line(colour = "black", 
                                   linetype = "solid")) +
  
  theme(legend.position="top" , legend.title=element_blank()) +
  
  theme(text=element_text(size=16, face = "bold")) +
  
  scale_color_manual(values=c("#5D6D7E","#D11300" )) + 
  
  theme(  panel.background = element_rect(fill = "white",
                                          colour = "white")) + 
  scale_fill_discrete(name = "") +
  geom_hline(yintercept= 0.05, size=1.5, linetype="dashed", alpha=0.7) +
 
  xlab("Metabolites") +
  ylab("P_value_squared") + scale_y_sqrt() #for maximize scale of the data



ggsave("lm.jpeg",x, dpi = 300, width = 8, height = 7)


significant_met <- significants[which(p_values< 0.05),]

annotation <- plsma_annotation[which(rownames(significant_met) %in% plsma_annotation$BIOCHEMICAL ), ]
  
annotation$BIOCHEMICAL



#--------------------------------------------------------------------------------------
