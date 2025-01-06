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
cor_xy_fla <- read_xlsx("../data/cor/GC-MC-Cor.xlsx")
names1 <- cor_xy_fla$Name
cor_xy_fla <- cor_xy_fla[-1]
rownames(cor_xy_fla) <- names1

cor_xy <- read_tsv("../data/cor/taxa_summary.tsv")
names2 <- cor_xy$genus
cor_xy <- cor_xy[-1]
rownames(cor_xy) <- names2

cor_xy_yx <- read.csv("../data/cor/yexiang_data.csv", row.names = 1)

cor_xy_1 <- apply(cor_xy, 2, as.numeric) %>% as.data.frame()#R语言将字符串矩阵转化为数值型矩阵 
rownames(cor_xy_1) <- rownames(cor_xy)

cor_xy_yx_1 <- apply(cor_xy_yx, 2, as.numeric) %>% as.data.frame()#R语言将字符串矩阵转化为数值型矩阵 
rownames(cor_xy_yx_1) <- rownames(cor_xy_yx)

cor_xy_fla_1 <- apply(cor_xy_fla, 2, as.numeric) %>% as.data.frame()#R语言将字符串矩阵转化为数值型矩阵 
rownames(cor_xy_fla_1) <- rownames(cor_xy_fla)

# cor.result<-corr.test(cor_xy_1,cor_xy_fla_1,method = "spearman")
# cor.result$r
# col3 <- grDevices::colorRampPalette(c("#2082DD", "white", "#FF3F3F")) 
# col3
# 用于绘制因素之间的相关性，可以将多种绘制方式混合运用，上下采用不同的展示方法，见于NC上的文章。
# corrplot.mixed(corr=corrmatrix,
#                lower="number",
#                upper="circle",
#                diag="u",
#                upper.col =col3(20),
#                lower.col = col3(20),
#                number.cex=0.9,
#                p.mat= res1$p,
#                sig.level= 0.05,
#                bg = "white",
#                is.corr = TRUE, 
#                outline = FALSE, 
#                mar = c(0,0,3,0),
#                addCoef.col = NULL, 
#                addCoefasPercent = FALSE, 
#                order = c("original"),
#                rect.col = "black", 
#                rect.lwd = 1, 
#                tl.cex = 1.2,
#                tl.col = "black", 
#                tl.offset = 0.4, 
#                tl.srt = 90,
#                cl.cex = 1.1, 
#                cl.ratio = 0.2,
#                tl.pos="lt",
#                cl.offset = 0.5 )
# res1 <-corrplot::cor.mtest(corrmatrix, conf.level= .95)
# res1$p
# res1$lowCI
# res1$uppCI

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
