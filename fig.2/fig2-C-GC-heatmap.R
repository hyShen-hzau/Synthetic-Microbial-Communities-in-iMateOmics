library(ComplexHeatmap)
library(readr)
library(tidyverse)
library(vegan)
library(circlize)
library(readxl)

data1 <- read_xlsx("sample-6type.xlsx")
group <- read_xlsx("Group.xlsx")
data1 <- data1[3:21]
group <- group[3:4]
colnames(data1)[1] <- "id"
colnames(group)[1] <- "id"
data1 <- data1 %>% column_to_rownames("id")
colSums(data1[1:6,])
x <- log2(data1+1)
#x <- decostand(data1,method="log") 
colSums(x[1:6,])
library(circlize)
max(x)
min(x)
col_fun <-  colorRamp2(c(0,1,2,3,4, max(x)), c("white","#ffebcd", "#80D1C8","#01b3bf", "#FF770F","#fe4c40"))
col_fun(seq(0, max(x)))

group1 <- data_frame(id=colnames(data1),
                     TYPE=c(rep("CK",3),rep("LAB",3),rep("WESS",3),rep("PICHI",3),rep("CAN",3),rep("SAC",3)))
group1$TYPE <- factor(group1$TYPE,levels = c("CK","LAB","WESS","PICHI","CAN" ,"SAC"),
                      ordered = T
)



p1 <- Heatmap(x,cluster_columns = F,cluster_rows = F,
              name = " ",col = col_fun,
              column_title_rot = 0,row_title_rot = 90,#设置标题旋转角度 
              row_names_gp = gpar(fontsize = 10), 
              column_names_gp = gpar(fontsize = 8),
              rect_gp = gpar(col = "white", lty = 1 ,lwd = 1),
              column_names_centered = F,
              na_col = "black", #若矩阵中存在NA值，特定设置为黑色
              border_gp = gpar(col = "black", lwd = 2), #设置热图的边框与粗度
              width = ncol(x)*unit(4.5, "mm"), 
              height = nrow(x)*unit(5, "mm"),#控制每一个单元格的大小
              row_title_gp = gpar(fontsize = 10, # 行标题大小
                                  ineheight = 1.2,
                                  #fontface = "bold" , # 行标题字体
                                  fill = c("#ffa39e","#ffd591","#fffb8f","#91caff","#ffadd2"),  # 行标题背景色
                                  col = "black", # 行标题颜色
                                  border = "black"# 行标题边框色
              )  ,
              row_split = group$Type   ,   #按照Type对热图行进行分割
              column_title_gp = gpar(fontsize = 8, # 行标题大小
                                     ineheight = 1,
                                     #fontface = "bold" , # 行标题字体
                                     fill = c("#ffffb8","#fffb8f","#fff566","#ffec3d","#d4b106","#ad8b00"),  # 行标题背景色
                                     col = "black", # 行标题颜色
                                     border = "black"# 行标题边框色
              )  ,
              column_split = group1$TYPE ,
              heatmap_legend_param = list(title = "log(x+1)", at = c(0,2,4,6,8))
              #row_km = 2按聚类行分隔
              #column 按列分隔
              
)
p1

pdf(file="result/sample-62.pdf")
p1
dev.off()

png(file="result/sample-611.png")
p1
dev.off()



#ck-lp----------------------------------------------------------
data2 <- read_xlsx("lp-ck-SANMPLE.xlsx")
group <- read_xlsx("Group.xlsx")
data2 <- data2[3:33]
group <- group[3:4]
colnames(data2)[1] <- "id"
colnames(group)[1] <- "id"
data2 <- data2 %>% column_to_rownames("id")

y <- log2(data2+1)
#y <- decostand(data2,method="hellinger") 
colSums(y)
rowSums(y)/30

library(circlize)
max(y)
min(y)
col_fun1 <-  colorRamp2(c(0,0.5,1,2,3,4,5,6, max(x)), 
                       c("#fff0f6","#ffd6e7","#ffadd2",
                         "#bae0ff","#69b1ff","#1677ff","#003eb3",
                         "#ffa940","#d46b08"))
col_fun1(seq(0, max(x)))

group2 <- data_frame(id=colnames(data2),
                        TYPE=c(rep("0 day",6),rep("3 day",6),rep("10 day",6),rep("30 day",6),rep("45 day",6)))
group2$TYPE <- factor(group2$TYPE,levels = c("0 day","3 day","10 day","30 day","45 day"),
                      ordered = T
                      )

p2 <- Heatmap(y,cluster_columns = F,cluster_rows = F,
              name = " ",col = col_fun1,
              column_title_rot = 0,row_title_rot = 90,#设置标题旋转角度 
              row_names_gp = gpar(fontsize = 10), 
              column_names_gp = gpar(fontsize = 8),
              rect_gp = gpar(col = "white", lty = 1 ,lwd = 1),
              column_names_centered = F,
              na_col = "black", #若矩阵中存在NA值，特定设置为黑色
              border_gp = gpar(col = "black", lwd = 2), #设置热图的边框与粗度
              width = ncol(x)*unit(4.5, "mm"), 
              height = nrow(x)*unit(5, "mm"),#控制每一个单元格的大小
              row_title_gp = gpar(fontsize = 10, # 行标题大小
                                  ineheight = 1.2,
                                  #fontface = "bold" , # 行标题字体
                                  fill = c("#ffa39e","#ffd591","#fffb8f","#91caff","#ffadd2"),  # 行标题背景色
                                  col = "black", # 行标题颜色
                                  border = "black"# 行标题边框色
              )  ,
              row_split = group$Type   ,   #按照Type对热图行进行分割
              column_title_gp = gpar(fontsize = 10, # 行标题大小
                                  ineheight = 1.5,
                                  #fontface = "bold" , # 行标题字体
                                  fill = c("#f6ffed","#d9f7be","#95de64","#52c41a","#237804"),  # 行标题背景色
                                  col = "black", # 行标题颜色
                                  border = "black"# 行标题边框色
              )  ,
              column_split = group2$TYPE ,
              heatmap_legend_param = list(title = "log(x+1)", at = c(0,0.5,1,2,3,4,5,6,7,8))
              #row_km = 2按聚类行分隔
              #column 按列分隔
              
)
p2

pdf(file="result/sample-lp-ck.pdf",width = 8.5, height = 9.05)
p2
dev.off()
pdf(file="result/sample-611.tiff")
p1
dev.off()