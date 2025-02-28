library(linkET)
library(ggplot2)
library(dplyr)
library(corrplot)

#读取示例数据
micro <- read.csv('xy.csv', row.names = 1)
env <- read.csv('fla_sta.csv', row.names = 1)

title <- colnames(micro)
micro <- t(micro) %>% as.data.frame()
rownames(env) <- title


#计算 Mantel 相关性
#通过 spec_select 指定数据组范围，例如这个示例数据中，微生物矩阵的第1-22列是物种丰度数据（指定名称 Taxonomy），第23-40列是基因丰度数据（指定名称 Function）
#默认情况下，对 spec 计算 Bray-Curtis 距离，对 env 计算 Euclidean 距离，然后计算二者 Mantel 相关
mantel <- mantel_test(spec = micro, env = env, 
                      spec_select = list(Bacteria = 1:6, Fungus = 7:10), 
                      mantel_fun = 'mantel')

mantel <- mutate(mantel, 
                 rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf), labels = c('< 0.2', '0.2 - 0.4', '>= 0.4')),
                 pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c('< 0.01', '0.01 - 0.05', '>= 0.05'))
)
mantel




qcorrplot(correlate(env, method = 'spearman'), type = 'upper', diag = FALSE) +  
  geom_square() +  #绘制 Spearman 相关系数热图
  #geom_mark(sep = '\n', size = 2.5, sig.thres = 0.05) +  #显示 Spearman 相关系数和显著性
  geom_couple(aes(color = pd, size = rd), data = mantel, curvature = nice_curvature()) +  #环境和微生物的相关性展示为上述 Mantel 相关
  scale_fill_gradientn(colors = c('#053061', '#68A8CF', 'white', '#F7B394', '#67001F'), limits = c(-1, 1)) + 
  scale_size_manual(values = c(0.5, 1, 2)) +  
  scale_color_manual(values = c('#D95F02', '#1B9E77', '#E0E0E0')) +  
  guides(color = guide_legend(title = "Mantel's p", order = 1), #图例标题和排序
         size = guide_legend(title = "Mantel's r", order = 2), 
         fill = guide_colorbar(title = "Spearman's r", order = 3)) +
  theme(legend.key = element_blank())



q <- qcorrplot(correlate(env, method = 'spearman'), type = 'upper', diag = FALSE) +  
  geom_square(percent = TRUE , rfill = "#F2AD00")+  #绘制 Spearman 相关系数热图
  #geom_mark(sep = '\n', size = 3, sig.thres = 0.05) +  #显示 Spearman 相关系数和显著性
  geom_couple(aes(color = pd, size = rd), data = mantel, curvature = nice_curvature()) +  #环境和微生物的相关性展示为上述 Mantel 相关
  scale_fill_gradientn(colors = c('#053061', '#68A8CF', 'white', '#F7B394', '#67001F'), limits = c(-1, 1)) + 
  scale_size_manual(values = c(0.5, 1.5, 2.5)) +  
  scale_color_manual(values = c('#D95F02', '#1B9E77', '#E0E0E0')) +  
  guides(color = guide_legend(title = "Mantel's p", order = 1), #图例标题和排序
         size = guide_legend(title = "Mantel's r", order = 2), 
         fill = guide_colorbar(title = "Spearman's r", order = 3)) +
  theme(legend.key = element_blank())


env <- data.frame(env)
result.cor <- cor(env, method = 'spearman')

ggsave(q,device = "pdf",file="result/mental_2024.pdf",width = 21, height = 18, units = "cm",onefile=F)  
ggsave(q,device = "png",file="result/mental_2024.png",width = 21, height = 18, units = "cm")
