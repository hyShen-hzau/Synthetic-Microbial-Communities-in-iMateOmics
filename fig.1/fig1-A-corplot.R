library(tidyverse)
library(permute)
library(lattice)
library(vegan)
library(phyloseq)
library(tidyverse)
library(pheatmap)
library(corrplot)


cor_xy_fla <- read.table("result/newcorplot/cor_xy_fla.txt",row.names =  1,sep = '\t')
colnames(cor_xy_fla) <- cor_xy_fla[1,]
cor_xy_fla <- cor_xy_fla[-1,]


cor_xy <- read.table("result/newcorplot/cor_xy_xy.txt",row.names =  1,sep = '\t')
colnames(cor_xy) <- cor_xy[1,]
cor_xy <- cor_xy[-1,]

cor_xy_ami <- read.table("result/newcorplot/cor_xy_ami.txt",row.names =  1,sep = '\t')
colnames(cor_xy_ami) <- cor_xy_ami[1,]
cor_xy_ami <- cor_xy_ami[-1,]


cor_xy_1 <- apply(cor_xy, 2, as.numeric) #R语言将字符串矩阵转化为数值型矩阵 
rownames(cor_xy_1) <- rownames(cor_xy)
corrplot(cor_xy_1, method = 'number', number.cex = 0.8, diag = FALSE, tl.cex = 0.8)
corrplot(cor_xy_1, add = TRUE, type = 'upper', method = 'pie', diag = FALSE, tl.pos = 'n', cl.pos = 'n',tl.col = "black")

corrplot(cor_xy_1, method = "pie", 
         type = "upper",
         tl.col = "black", tl.cex = 0.9, tl.srt = 45
)



cor_xy_ami_1 <- apply(cor_xy_ami, 2, as.numeric) #R语言将字符串矩阵转化为数值型矩阵 
rownames(cor_xy_ami_1) <- rownames(cor_xy_ami)


corrplot(cor_xy_ami_1, method = "pie", 
         type = "full",
         tl.col = "black", tl.cex = 0.6, tl.srt = 45
)


cor_xy_fla_1 <- apply(cor_xy_fla, 2, as.numeric) #R语言将字符串矩阵转化为数值型矩阵 
rownames(cor_xy_fla_1) <- rownames(cor_xy_fla)


corrplot(cor_xy_fla_1, method = "pie", 
         type = "full",
         tl.col = "black", tl.cex = 0.8, tl.srt = 45
)



#test----------------------------------------------
library(tidyverse)
library(readxl)
library(magrittr)
df2 <- cor_xy_ami_1 %>% as.data.frame() %>% rownames_to_column(var="ID") %>%  
  pivot_longer(-1) %>% set_colnames(c("ID","name","value"))
df2$ID <- factor(df2$ID,levels = rownames(cor_xy_ami),ordered = T)

df3 <- cor_xy_fla_1 %>% as.data.frame() %>% rownames_to_column(var="ID") %>%  
  pivot_longer(-1) %>% set_colnames(c("ID","name","value"))
df3$ID <- factor(df3$ID,levels = rownames(cor_xy_fla),ordered = T)
# p1 <- df2 %>% filter(name %in% c("J043V6","J035V8","J009V6")) %>% 
#   mutate(p_value=case_when(value <0 ~ "**")) %>% drop_na() %>% 
#   mutate(group=case_when(value < -0.5 ~ "latter",
#                          TRUE ~ "former"))

df4 <- rbind(df2,df3)
x1 <- unique(df4$ID)
df4$ID <- factor(df4$ID,levels = x1,ordered = T)



ggplot()+ 
  geom_tile(data=df4  ,
            aes(name,ID,fill=value))+
  # geom_text(data=df2,aes(name, ID,label=p_value),
  #           color="white",vjust=1,size=5,hjust=0.5)+ 
  geom_point(data=df4 %>% 
               filter(ID %in% rownames(cor_xy_fla)) ,
             aes(name,factor(ID,levels = x1),color=value,size=value))+
  coord_cartesian(clip = "off") +  # 关闭坐标轴裁剪
  labs(x=NULL,y=NULL)+  # 移除x和y轴标签
  scale_x_discrete(expand=c(0,0))+  # 设置x轴
  scale_fill_gradient2(low="#788FCE",mid="#FBFEF9",high="#E6956F")+
  scale_color_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446")+  # 设置填充颜色
  theme_test()+  # 使用theme_test主题
  theme(axis.text.x=element_text(color="black",size=8,face="bold",angle = 45, hjust = 1),  # 自定义x轴文本
        axis.text.y=element_text(color="black",size=8,face="bold",angle = 0,vjust=0.5),  # 自定义y轴文本
        axis.ticks = element_blank(),  # 移除轴刻度
        legend.title = element_blank(),  # 移除图例标题
        legend.background = element_blank(),  # 设置图例背景为透明
        legend.text = element_text(size=8,color="black",face="bold"),  # 自定义图例文本
        legend.position ="right",  # 设置图例位置
        legend.spacing.x = unit(0.05,"in"),  # 设置图例内部的水平间距
        plot.margin = ggplot2::margin(10,30,10,10),
        legend.direction = "vertical",
        legend.box = "horizontal")




ggplot()+ 
  geom_tile(data=df2 ,
             aes(name,ID,fill=value))+
  geom_point(data=df4 %>% 
              filter(ID %in% rownames(cor_xy_fla)) ,
             aes(name,factor(ID,levels = x1),color=value,size=value))+
  coord_cartesian(clip = "off") +  # 关闭坐标轴裁剪
  labs(x=NULL,y=NULL)+  # 移除x和y轴标签
  scale_x_discrete(expand=c(0,0))+  # 设置x轴
  scale_fill_gradient2(low="#788FCE",mid="#FBFEF9",high="#E6956F")+
  scale_color_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446")+  # 设置填充颜色
  theme_test()+  # 使用theme_test主题
  theme(axis.text.x=element_text(color="black",size=8,face="bold",angle = 45, hjust = 1),  # 自定义x轴文本
        axis.text.y=element_text(color="black",size=8,face="bold",angle = 0,vjust=0.5),  # 自定义y轴文本
        axis.ticks = element_blank(),  # 移除轴刻度
        legend.title = element_blank(),  # 移除图例标题
        legend.background = element_blank(),  # 设置图例背景为透明
        legend.text = element_text(size=8,color="black",face="bold"),  # 自定义图例文本
        legend.position ="right",  # 设置图例位置
        legend.spacing.x = unit(0.05,"in"),  # 设置图例内部的水平间距
        plot.margin = ggplot2::margin(10,30,10,10),
        legend.direction = "vertical",
        legend.box = "horizontal")


#ggplot------------------------------------------------------
library(tidyverse)
library(readxl)
library(magrittr)

rownames(cor_xy_ami_1) <- str_c(".",rownames(cor_xy_ami_1))
rownames(cor_xy_fla_1)[11] <- "alpha-Terpineol"
dd4 <- rbind(cor_xy_ami_1,cor_xy_fla_1)
df41 <- dd4 %>% as.data.frame() %>% rownames_to_column(var="ID") %>%  
  pivot_longer(-1) %>% set_colnames(c("ID","name","value"))
df41$ID <- factor(df41$ID,levels = rownames(dd4),ordered = T)
p1 <- df41 %>% filter(ID %in% rownames(cor_xy_ami_1)) 

q1 <- 
  ggplot()+ 
  geom_tile(data=p1 ,
            aes(name,ID,fill=value))+
  # geom_text(data=df2,aes(name, ID,label=p_value),
  #           color="white",vjust=1,size=5,hjust=0.5)+ 
  geom_point(data=df41 %>% 
               filter(ID %in% rownames(cor_xy_fla)) ,
             aes(name,ID,color=value,size=value))+
  coord_cartesian(clip = "off") +  # 关闭坐标轴裁剪
  labs(x=NULL,y=NULL)+  # 移除x和y轴标签
  scale_x_discrete(expand=c(0,0))+  # 设置x轴
  scale_fill_gradient2(low="#788FCE",mid="#FBFEF9",high="#E6956F")+
  scale_color_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446")+  # 设置填充颜色
  theme_test()+  # 使用theme_test主题
  theme(axis.text.x=element_text(color="black",size=8,face="bold",angle = 45, hjust = 1),  # 自定义x轴文本
        axis.text.y=element_text(color="black",size=8,face="bold",angle = 0,vjust=0.5),  # 自定义y轴文本
        axis.ticks = element_blank(),  # 移除轴刻度
        legend.title = element_blank(),  # 移除图例标题
        legend.background = element_blank(),  # 设置图例背景为透明
        legend.text = element_text(size=8,color="black",face="bold"),  # 自定义图例文本
        legend.position ="right",  # 设置图例位置
        legend.spacing.x = unit(0.05,"in"),  # 设置图例内部的水平间距
        plot.margin = ggplot2::margin(10,30,10,10),
        legend.direction = "vertical",
        legend.box = "horizontal")
ggsave(q1,device = "pdf",file="result/cor-2023/cor_tol.pdf",width = 21, height = 16, units = "cm",onefile=F)   
ggsave(q1,device = "png",file="result/cor-2023/cor_tol.png",
       dpi = 300) 


p1 <- ggplot()+ 
  geom_tile(data=df2 ,
            aes(name,ID,fill=value))+
  labs(x=NULL,y=NULL)+  # 移除x和y轴标签
  scale_x_discrete(expand=c(0,0))+  # 设置x轴
  scale_fill_gradient2(low="#788FCE",mid="#FBFEF9",high="#E6956F")+
  theme_test()+  # 使用theme_test主题
  theme(axis.text.x=element_blank(),  # 自定义x轴文本
        axis.text.y=element_text(color="black",size=8,face="bold",angle = 0,vjust=0.5),  # 自定义y轴文本
        axis.ticks = element_blank(),  # 移除轴刻度
        legend.title = element_blank(),  # 移除图例标题
        legend.background = element_blank(),  # 设置图例背景为透明
        legend.text = element_text(size=8,color="black",face="bold"),  # 自定义图例文本
        legend.position ="right",  # 设置图例位置
        legend.spacing.x = unit(0.05,"in"),  # 设置图例内部的水平间距
        plot.margin = ggplot2::margin(10,30,10,10),
        legend.direction = "vertical",
        legend.box = "horizontal")
ggsave(p1,device = "pdf",file="result/cor-2023/cor_ami.pdf",width = 21, height = 16, units = "cm",onefile=F)   
ggsave(p1,device = "png",file="result/cor-2023/cor_ami.png",
       dpi = 300) 


p2 <- ggplot()+ 
  geom_point(data=df4 %>% 
               filter(ID %in% rownames(cor_xy_fla)) ,
             aes(name,ID,color=value,size=value))+
  coord_cartesian(clip = "off") +  # 关闭坐标轴裁剪
  labs(x=NULL,y=NULL)+  # 移除x和y轴标签
  scale_x_discrete(expand=c(0,0))+  # 设置x轴
  scale_color_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446")+  # 设置填充颜色
  theme_test()+  # 使用theme_test主题
  theme(axis.text.x=element_text(color="black",size=8,face="bold",angle = 45, hjust = 1),  # 自定义x轴文本
        axis.text.y=element_text(color="black",size=8,face="bold",angle = 0,vjust=0.5),  # 自定义y轴文本
        axis.ticks = element_blank(),  # 移除轴刻度
        legend.title = element_blank(),  # 移除图例标题
        legend.background = element_blank(),  # 设置图例背景为透明
        legend.text = element_text(size=8,color="black",face="bold"),  # 自定义图例文本
        legend.position ="right",  # 设置图例位置
        legend.spacing.x = unit(0.05,"in"),  # 设置图例内部的水平间距
        plot.margin = ggplot2::margin(10,30,10,10),
        legend.direction = "vertical",
        legend.box = "horizontal")

ggsave(p2,device = "pdf",file="result/cor-2023/cor_fla.pdf",width = 21, height = 16, units = "cm",onefile=F)   
ggsave(p2,device = "png",file="result/cor-2023/cor_fla.png",
       dpi = 300) 
