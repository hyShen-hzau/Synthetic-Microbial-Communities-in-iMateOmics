library(readr)
library(tidyverse)
#首先加载上述自定义函数 raup_crick()，不多说
source('RaupCrick.txt')

#然后读取示例的分类群丰度表
#上述自定义函数 raup_crick() 要求输入矩阵里面行是样本，列是物种或分类群
spp <- read.csv('../data/bac-fungal/bac.csv', row.names = 1) %>% column_to_rownames("taxonomy") %>% t() %>% 
  as.data.frame() %>% rownames_to_column("X") %>% remove_rownames()

#使用上述自定义函数 raup_crick() 计算所有样本（群落）对之间的 Raup-Crick 相异指数
#详细参数设置，请审阅源代码里面的注释
set.seed(123)
raup_crick.dist <- raup_crick(spXsite = spp, 
    plot_names_in_col1 = TRUE, 
    classic_metric = FALSE, 
    split_ties = TRUE, 
    reps = 999, 
    set_all_species_equal = FALSE, 
    as.distance.matrix = TRUE, 
    report_similarity = FALSE
)

#计算好的 Raup-Crick 相异指数默认以 dist 类型存储，可转换为 matrix 类型输出
raup_crick.matrix <- as.matrix(raup_crick.dist)
#write.table(raup_crick.matrix, 'Raup-Crick-bac.txt', sep = '\t', col.names = NA)

# raup_crick.matrix <- read.table('../result/beta_rc/Raup-Crick-bac.txt', row.names = 1)
# colnames(raup_crick.matrix) <- raup_crick.matrix[1,]
# raup_crick.matrix <- raup_crick.matrix[-1,]

##基于 Raup-Crick 指数的 NMDS 分析的例子
library(vegan)
library(ggplot2)

#读取计算好的 Raup-Crick 相异指数矩阵
beta_RC <- read.delim('../data/Raup-Crick-bac.txt', row.names = 1, sep = '\t')
beta_RC.dist <- as.dist(beta_RC)	#转为 dist 数据类型

#NMDS 排序，定义 2 个维度，详情 ?metaMDS
nmds_dis <- metaMDS(beta_RC.dist, k = 2)
nmds_dis

#获取 stress 值
stress <- nmds_dis$stress
stress

#提取样本得分（排序坐标）
nmds_dis_site <- data.frame(nmds_dis$points)
nmds_dis_site

#添加分组信息，绘制 NMDS 图
nmds_dis_site$group <- c(rep('CK', 15), rep('Lpscw', 15))
nmds_dis_site$time <- c(rep('0day', 3), rep('3day', 3),rep('10day', 3), rep('30day', 3),rep('45day', 3),
                        rep('0day', 3), rep('3day', 3),rep('10day', 3), rep('30day', 3),rep('45day', 3))

p1 <- ggplot(data = nmds_dis_site, aes(MDS1, MDS2)) +
geom_point(aes(color = group, shape = time),size = 8) +
scale_shape_manual(values = c(0,1,2,5,11)) +
scale_color_manual(values = c('red3', 'orange3')) +
scale_fill_manual(values = c('red', 'orange')) +
stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) +	
theme(panel.grid.major = element_line(color = 'gray', size = 0.2), panel.background = element_rect(color = 'black', fill = 'transparent')) +
labs(x = 'NMDS1', y = 'NMDS2') +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = rel(1.5),face = "italic"),
        legend.text = element_text(size = 16),
        plot.title = element_text(
          size = 24,        # 标题字体大小
          face = "bold",    # 标题字体样式 ("plain"普通, "bold"加粗, "italic"斜体, "bold.italic"加粗斜体)
          color = "black",   # 标题字体颜色
          hjust = 0.5       # 标题水平对齐 (0为左对齐，0.5为居中，1为右对齐)
        )
  )+
annotate('text', label = paste('Stress =', round(stress, 5)), x = 0.64, y = 0.6, size = 6, colour = 'black') +
  ggtitle("Bacterial community")
p1


env1 <- as.dist(beta_RC[1:15,1:15])
env2 <- as.dist(beta_RC[16:30,16:30])

plot_data <- data.frame(
    beta_RC = c(env1, env2), 
    group = c(rep('CK', length(env1)), rep('Lpscw', length(env2)))
)

p2 <- ggplot(plot_data, aes(group, beta_RC, fill = group)) +
geom_boxplot(show.legend = FALSE) +
scale_fill_manual(values = c('red', 'orange')) +
labs(x = '', y = 'Raup-Crick Dissimilarity') +
  annotate('text', label = '45d2:3d1-', x = 0.89, y = 0.82, size = 4, colour = 'red') +
  annotate('text', label = '45d2:3d2-', x = 0.89, y = 0.96, size = 4, colour = 'orange')+
  theme_bw()+
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = rel(1.5),face = "italic"),
        legend.text = element_text(size = 16),
        plot.title = element_text(
          size = 24,        # 标题字体大小
          face = "bold",    # 标题字体样式 ("plain"普通, "bold"加粗, "italic"斜体, "bold.italic"加粗斜体)
          color = "black",   # 标题字体颜色
          hjust = 0.5       # 标题水平对齐 (0为左对齐，0.5为居中，1为右对齐)
        )
  )+
  ggtitle("Bacterial community")
p2




#its----------------------------------------------------------
##基于 Raup-Crick 指数的 NMDS 分析的例子
library(vegan)
library(ggplot2)

#读取计算好的 Raup-Crick 相异指数矩阵
beta_RC <-  read.delim('../data/Raup-Crick-its.txt', row.names = 1, sep = '\t')
beta_RC.dist <- as.dist(beta_RC)	#转为 dist 数据类型

#NMDS 排序，定义 2 个维度，详情 ?metaMDS

nmds_dis <- metaMDS(beta_RC.dist, k = 2)
nmds_dis

#获取 stress 值
stress <- nmds_dis$stress
stress

#提取样本得分（排序坐标）
nmds_dis_site <- data.frame(nmds_dis$points)
nmds_dis_site

#添加分组信息，绘制 NMDS 图
nmds_dis_site$group <- c(rep('CK', 15), rep('Lpscw', 15))
nmds_dis_site$time <- c(rep('0day', 3), rep('3day', 3),rep('10day', 3), rep('30day', 3),rep('45day', 3),
                        rep('0day', 3), rep('3day', 3),rep('10day', 3), rep('30day', 3),rep('45day', 3))

p3 <- ggplot(data = nmds_dis_site, aes(MDS1, MDS2)) +
  geom_point(aes(color = group, shape = time),size = 8) +
  scale_shape_manual(values = c(0,1,2,5,11)) +
  scale_color_manual(values = c('red3', 'orange3')) +
  scale_fill_manual(values = c('red', 'orange')) +
  stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) +	
  theme(panel.grid.major = element_line(color = 'gray', size = 0.2), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  labs(x = 'NMDS1', y = 'NMDS2')+
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = rel(1.5),face = "italic"),
        legend.text = element_text(size = 16),
        plot.title = element_text(
          size = 24,        # 标题字体大小
          face = "bold",    # 标题字体样式 ("plain"普通, "bold"加粗, "italic"斜体, "bold.italic"加粗斜体)
          color = "black",   # 标题字体颜色
          hjust = 0.5       # 标题水平对齐 (0为左对齐，0.5为居中，1为右对齐)
        )
  )+
  ggtitle("Fungal community")+
  annotate('text', label = paste('Stress =', round(stress, 3)), x = 0.35, y = 0.4, size = 6, colour = 'black') 
p3



##再象征性画个 Raup-Crick 相异指数的组内箱线图
env1 <- as.dist(beta_RC[1:15,1:15])
env2 <- as.dist(beta_RC[16:30,16:30])

plot_data <- data.frame(
  beta_RC = c(env1, env2), 
  group = c(rep('CK', length(env1)), rep('Lpscw', length(env2)))
)

p4 <- ggplot(plot_data, aes(group, beta_RC, fill = group)) +
  geom_boxplot(show.legend = FALSE) +
  scale_fill_manual(values = c('red', 'orange')) +
  labs(x = '', y = 'Raup-Crick Dissimilarity') +
  annotate('text', label = '45d1:0d3-', x = 1.89, y = 0.94, size = 4, colour = 'red') +
  annotate('text', label = '-45d2:0d3', x = 2.11, y = 0.94, size = 4, colour = 'orange')+
  annotate('text', label = '45d3:0d3-', x = 1.89, y = 0.79, size = 4, colour = 'black')+
  theme_bw()+
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = rel(1.5),face = "italic"),
        legend.text = element_text(size = 16),
        plot.title = element_text(
          size = 24,        # 标题字体大小
          face = "bold",    # 标题字体样式 ("plain"普通, "bold"加粗, "italic"斜体, "bold.italic"加粗斜体)
          color = "black",   # 标题字体颜色
          hjust = 0.5       # 标题水平对齐 (0为左对齐，0.5为居中，1为右对齐)
        )
  )+
  ggtitle("Fungal community")
p4
ggsave(p4, device = "pdf", filename = "../result/beta_rc/its_box-1_rev.pdf",
       width = 15, height = 10, units = "cm")
ggsave(p4, device = "jpg", filename = "../result/beta_rc/its_box-1_rev.jpg",
       width = 9, height = 6, units = "in", dpi = 300)
 