library(tidyverse)
library(permute)
library(lattice)
library(vegan)
library(phyloseq)
library(tidyverse)
library(pheatmap)
library(corrplot)
library(readxl)
library(psych)


#数据导入------------------
cor_xy_fla <- read_xlsx("cor/GC-MC-Cor.xlsx")
names1 <- cor_xy_fla$Name
cor_xy_fla <- cor_xy_fla[-1]
rownames(cor_xy_fla) <- names1

cor_xy <- read_tsv("cor/taxa_summary.tsv")
names2 <- cor_xy$genus
cor_xy <- cor_xy[-1]
rownames(cor_xy) <- names2

cor_xy_yx <- read.csv("cor/yexiang_data.csv", row.names = 1)

cor_xy_1 <- apply(cor_xy, 2, as.numeric) %>% as.data.frame()#R语言将字符串矩阵转化为数值型矩阵 
rownames(cor_xy_1) <- rownames(cor_xy)

cor_xy_yx_1 <- apply(cor_xy_yx, 2, as.numeric) %>% as.data.frame()#R语言将字符串矩阵转化为数值型矩阵 
rownames(cor_xy_yx_1) <- rownames(cor_xy_yx)

cor_xy_fla_1 <- apply(cor_xy_fla, 2, as.numeric) %>% as.data.frame()#R语言将字符串矩阵转化为数值型矩阵 
rownames(cor_xy_fla_1) <- rownames(cor_xy_fla)


library(RColorBrewer)#用RColorBrewer提供配色方案
coul = colorRampPalette(c("#55aeda","white","#fa6464"))(20)

#corrmatrix <- cor(t(cor_xy_1[1:15]),t(cor_xy_fla_1[1:15]),method = "spearman")

#CK-xy-fla-results
corr <- corr.test(t(cor_xy_1[1:15]),t(cor_xy_fla_1[1:15]),method = "spearman")
corrmatrix <- corr$r

pdf(height=6.0, width=6.5,family='serif', bg = "white",
    file = "../result/corrplot/cor_fla_ck.pdf")#设置导出的图的高度、宽度、字体、背景色、格式
corrplot(corrmatrix,col= coul,order = "original",method = "ellipse",tl.cex= 0.9,tl.col ="black",
         tl.srt = 45, cl.align.text = "l", cl.length = 7, cl.ratio = 0.3, cl.cex = 0.9, 
         p.mat = corr$p, insig="label_sig", sig.level = c(0.001,0.01,0.05),
         pch.cex = 0.5,mar = c(0.5,0.2,0.2,0.2))
dev.off()


png(file='../result/corrplot/cor_fla_ck.png',family='serif', height=1950, 
    width=1500,units = "px",bg = "white",res = 300)
corrplot(corrmatrix,col= coul,order = "original",method = "ellipse",tl.cex= 0.9,tl.col ="black",
         tl.srt = 45, cl.align.text = "l", cl.length = 7, cl.ratio = 0.3, cl.cex = 0.9, 
         p.mat = corr$p, insig="label_sig", sig.level = c(0.001,0.01,0.05),
         pch.cex = 0.5,mar = c(0.5,0.2,0.2,0.2))
dev.off()


#Lp-xy-fla-results
corr1 <- corr.test(t(cor_xy_1[16:30]),t(cor_xy_fla_1[16:30]),method = "spearman")
corrmatrix1 <- corr1$r

pdf(height=6.0, width=6.5,family='serif', bg = "white",
    file = "../result/corrplot/cor_fla_lp.pdf")#设置导出的图的高度、宽度、字体、背景色、格式
corrplot(corrmatrix1,col= coul,order = "original",method = "ellipse",tl.cex= 0.9,tl.col ="black",
         tl.srt = 45, cl.align.text = "l", cl.length = 7, cl.ratio = 0.3, cl.cex = 0.9, 
         p.mat = corr1$p, insig="label_sig", sig.level = c(0.001,0.01,0.05),
         pch.cex = 0.5,mar = c(0.5,0.2,0.2,0.2))
dev.off()

png(file='../result/corrplot/cor_fla_lp.png',family='serif', height=1950, 
    width=1500,units = "px",bg = "white",res = 300)
corrplot(corrmatrix1,col= coul,order = "original",method = "ellipse",tl.cex= 0.9,tl.col ="black",
         tl.srt = 45, cl.align.text = "l", cl.length = 7, cl.ratio = 0.3, cl.cex = 0.9, 
         p.mat = corr1$p, insig="label_sig", sig.level = c(0.001,0.01,0.05),
         pch.cex = 0.5,mar = c(0.5,0.2,0.2,0.2))
dev.off()



#CK-xy-yx-results
corr2 <- corr.test(t(cor_xy_1[1:15]),t(cor_xy_yx_1[1:15]),method = "spearman")
corrmatrix2 <- corr2$r

pdf(height=6.0, width=6.5,family='serif', bg = "white",
    file = "../result/corrplot/cor_yx_ck.pdf")#设置导出的图的高度、宽度、字体、背景色、格式
corrplot(corrmatrix2,col= coul,order = "original",method = "ellipse",tl.cex= 0.9,tl.col ="black",
         tl.srt = 45, cl.align.text = "l", cl.length = 7, cl.ratio = 0.3, cl.cex = 0.9, 
         p.mat = corr2$p, insig="label_sig", sig.level = c(0.001,0.01,0.05),
         pch.cex = 0.5,mar = c(0.5,0.2,0.2,0.2))
dev.off()


png(file='../result/corrplot/cor_yx_ck.png',family='serif', height=1950, 
    width=1500,units = "px",bg = "white",res = 300)
corrplot(corrmatrix2,col= coul,order = "original",method = "ellipse",tl.cex= 0.9,tl.col ="black",
         tl.srt = 45, cl.align.text = "l", cl.length = 7, cl.ratio = 0.3, cl.cex = 0.9, 
         p.mat = corr2$p, insig="label_sig", sig.level = c(0.001,0.01,0.05),
         pch.cex = 0.5,mar = c(0.5,0.2,0.2,0.2))
dev.off()


#Lp-xy-yx-results
corr3 <- corr.test(t(cor_xy_1[16:30]),t(cor_xy_yx_1[16:30]),method = "spearman")
corrmatrix3 <- corr3$r

pdf(height=6.0, width=6.5,family='serif', bg = "white",
    file = "../result/corrplot/cor_yx_lp.pdf")#设置导出的图的高度、宽度、字体、背景色、格式
corrplot(corrmatrix3,col= coul,order = "original",method = "ellipse",tl.cex= 0.9,tl.col ="black",
         tl.srt = 45, cl.align.text = "l", cl.length = 7, cl.ratio = 0.3, cl.cex = 0.9, 
         p.mat = corr3$p, insig="label_sig", sig.level = c(0.001,0.01,0.05),
         pch.cex = 0.5,mar = c(0.5,0.2,0.2,0.2))
dev.off()


png(file='../result/corrplot/cor_yx_lp.png',family='serif', height=1950, 
    width=1500,units = "px",bg = "white",res = 300)
corrplot(corrmatrix3,col= coul,order = "original",method = "ellipse",tl.cex= 0.9,tl.col ="black",
         tl.srt = 45, cl.align.text = "l", cl.length = 7, cl.ratio = 0.3, cl.cex = 0.9, 
         p.mat = corr3$p, insig="label_sig", sig.level = c(0.001,0.01,0.05),
         pch.cex = 0.5,mar = c(0.5,0.2,0.2,0.2))
dev.off()
