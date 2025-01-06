##加载R包，初始传递命令
library(circlize) #使用该包绘制 circos 图
library(reshape2)
library(ComplexHeatmap) #可用此包添加图例
library(grid) #可用此包调整画板
library(tidyverse)
library(phyloseq)


#预处理--------------------------------------------------------------------
meta.kegg <- read_tsv("../data/function_compare/HG--YC--HS--JZ/meta.kegg_TPM.tsv")

name.kegg <- meta.kegg$Gene
name.kegg <- str_replace(name.kegg,"[|]",";")#x4
name.kegg <- str_replace(name.kegg,"[|]",";")
name.kegg <- str_replace(name.kegg,"[|]",";")
name.kegg <- str_replace(name.kegg,"[|]",";")
name.kegg
meta.kegg$taxonomy <- name.kegg
meta.kegg <- meta.kegg[-1]

rank1 <- map_dfr(meta.kegg$taxonomy, parse_taxonomy_qiime)$Rank1
rank2 <- map_dfr(meta.kegg$taxonomy, parse_taxonomy_qiime)$Rank2
rank3 <- map_dfr(meta.kegg$taxonomy, parse_taxonomy_qiime)$Rank3
rank4 <- map_dfr(meta.kegg$taxonomy, parse_taxonomy_qiime)$Rank4
unique(rank1)

text <- cbind(meta.kegg,rank1)
text <- cbind(text,rank4)
text <- text %>% filter(!is.na(rank4)) %>% select(-taxonomy)

r4 <- str_split(text$rank4,"[(]",n=2,simplify = T) %>% as.data.frame()
r4$V2 <- str_replace(r4$V2,"[)]","")
text$rank4 <- r4$V2

r1 <- str_split(text$rank1,"[(]",n=2,simplify = T) %>% as.data.frame()
text$rank1 <- r1$V1

text1 <- text %>% remove_rownames() %>% 
  group_by(rank1,rank4) %>%
  summarise(across(.cols = everything(),.fns = sum))


text2 <- text1  %>% 
  rowwise() %>% 
  mutate(total=sum(across(-1)),.after=1) %>% 
  ungroup() %>% 
  mutate(prop=total/sum(total),.after=1) %>% 
  arrange(desc(prop)) %>% 
  filter(row_number() <=20) %>% 
  select(-c("prop","total")) 


taxonomy <- text2[1:2]
colnames(taxonomy) <- c("OTU_ID",	"phylum")
unique(taxonomy$phylum)

write.csv(taxonomy,"../data/circlize/kegg/taxonomy.csv",row.names=FALSE)

group <- read.csv("../data/matedata/group_jz.csv") %>% as.data.frame()
topkegg <- text2 %>% select(-rank1) %>% t()
colnames(topkegg) <- topkegg[1,]
topkegg <- topkegg[-1,] %>% as.data.frame()
topkegg$Sampleid <- rownames(topkegg)

topkegg <- topkegg %>% remove_rownames()
topkegg$group <- group$Sampling_location


name <- topkegg$group
test <- topkegg[1:20]
test <- as.data.frame(lapply(topkegg, as.numeric))
test$group <- name

kg <- test  %>% group_by(group) %>% 
  summarise(across(.fns = mean))

topk <- t(kg) %>% as.data.frame()
colnames(topk) <- topk[1,]
topk <- topk[-1,]
topk <- topk[-21,]
write.csv(topk,"../data/circlize/kegg/topkegg.csv")







#kegg预处理---------------------------------------------------------------
##读取数据
#读取 taxonomy.txt 的内容，获取“OTU/分类”排序，OTU 名称
taxonomy <- read_csv("../data/circlize/kegg/taxonomy.csv")
tax_phylum <- unique(taxonomy$phylum)
taxonomy$phylum<- factor(taxonomy$phylum, levels = tax_phylum)
taxonomy <- taxonomy[order(taxonomy$phylum , taxonomy$OTU_ID), ]
all_otu <- taxonomy$OTU_ID

#读取 group.txt 的内容，获取“样本/分组”排序，样本名称
group <- read_csv("../data/circlize/kegg/group.csv")
all_group <- unique(group$group_ID)
group$group_ID <- factor(group$group_ID, levels = all_group)
all_sample <- group$sample_ID

#读取 otu_table.txt，排序 OTU 和样本
otu_table <- read_csv("../data/circlize/kegg/topkegg.csv")
otu_table <- merge(taxonomy, otu_table, by = 'OTU_ID')
otu_table <- otu_table[order(otu_table$phylum, otu_table$OTU_ID), ]
rownames(otu_table) <- otu_table$OTU_ID
otu_table <- otu_table[all_sample]

##生成作图数据
#circlize 外圈属性数据
all_ID <- c(all_otu, all_sample)
accum_otu <- rowSums(otu_table)
accum_sample <- colSums(otu_table)
all_ID_xlim <- cbind(rep(0, length(all_ID)),data.frame(c(accum_otu, accum_sample)))

#circlize 内圈连线数据
otu_table$otu_ID <- all_otu
plot_data <- melt(otu_table, id = 'otu_ID') #此处使用了reshape2包中的melt()命令
colnames(plot_data)[2] <- 'sample_ID'
plot_data$otu_ID <- factor(plot_data$otu_ID, levels = all_otu)
plot_data$sample_ID <- factor(plot_data$sample_ID, levels = all_sample)
plot_data <- plot_data[order(plot_data$otu_ID, plot_data$sample_ID), ]
plot_data <- plot_data[c(2, 1, 3, 3)]

#颜色设置
color_otu <- c('#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5', '#BC80BD', '#CCEBC5', '#FFED6F', '#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#66C2A5')
color_group <- c( "#E64B35E5", "#4DBBD5E5" ,"#00A087E5" ,"#3C5488E5")
color_phylum <- c('#BEAED4', '#FDC086', '#FFFF99', '#386CB0', '#F0027F',"#FBD26A","#01847f")

color_sample <- c( "#E64B35E5", "#4DBBD5E5" ,"#00A087E5" ,"#3C5488E5")

names(color_otu) <- all_otu
names(color_sample) <- all_sample

####circlize 绘图
pdf('circlize_plot.pdf', width = 20, height = 8)
circle_size = unit(1, 'snpc')

##整体布局
gap_size <- c(rep(3, length(all_otu) - 1), 6, rep(3, length(all_sample) - 1), 6)#定义间距
circos.par(cell.padding = c(0, 0, 0, 0), start.degree = 260, gap.degree = gap_size)#画板设置
circos.initialize(factors = factor(all_ID, levels = all_ID), xlim = all_ID_xlim)#定义绘图因子

##绘制 OTU 分类、样本分组区块（第一圈） 
circos.trackPlotRegion(
  ylim = c(0, 1), track.height = 0.05, bg.border = NA,
  panel.fun = function(x, y) {
    sector.index = get.cell.meta.data('sector.index')
    xlim = get.cell.meta.data('xlim')
    ylim = get.cell.meta.data('ylim')
  } )
 
for (i in 1:length(tax_phylum)) {
	tax_OTU <- {subset(taxonomy, phylum == tax_phylum[i])}$OTU_ID
	highlight.sector(tax_OTU, track.index = 1, col = color_phylum[i], text = tax_phylum[i], cex = 0.5, text.col = 'black', niceFacing = FALSE)
}

for (i in 1:length(all_group)) {
	group_sample <- {subset(group, group_ID == all_group[i])}$sample_ID
	highlight.sector(group_sample, track.index = 1, col = color_group[i], text = all_group[i], cex = 0.7, text.col = 'black', niceFacing = FALSE)
}

##各 OTU、样本绘制区
#添加百分比注释（第二圈）
circos.trackPlotRegion(
	ylim = c(0, 1), track.height = 0.05, bg.border = NA, 
	panel.fun = function(x, y) {
		sector.index = get.cell.meta.data('sector.index')
		xlim = get.cell.meta.data('xlim')
		ylim = get.cell.meta.data('ylim')
	} )

circos.track(
	track.index = 2, bg.border = NA, 
	panel.fun = function(x, y) {
		xlim = get.cell.meta.data('xlim')
		ylim = get.cell.meta.data('ylim')
		sector.name = get.cell.meta.data('sector.index')
		xplot = get.cell.meta.data('xplot')
		
		by = ifelse(abs(xplot[2] - xplot[1]) > 30, 0.25, 1)
		for (p in c(0, seq(by, 1, by = by))) circos.text(p*(xlim[2] - xlim[1]) + xlim[1], mean(ylim) + 0.4, paste0(p*100, '%'), cex = 0.4, adj = c(0.5, 0), niceFacing = FALSE)
		
		circos.lines(xlim, c(mean(ylim), mean(ylim)), lty = 3)
	} )

#绘制 OTU、样本主区块（第三圈）
circos.trackPlotRegion(
	ylim = c(0, 1), track.height = 0.03, bg.col = c(color_otu, color_sample), bg.border = NA, track.margin = c(0, 0.01),
	panel.fun = function(x, y) {
		xlim = get.cell.meta.data('xlim')
		sector.name = get.cell.meta.data('sector.index')
		circos.axis(h = 'top', labels.cex = 0.4, major.tick.percentage = 0.4, labels.niceFacing = FALSE)
		circos.text(mean(xlim), 0.2, sector.name, cex = 0.4, niceFacing = FALSE, adj = c(0.5, 0))
	} )

#绘制 OTU、样本副区块（第四圈）
circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.03, track.margin = c(0, 0.01))

##绘制 OTU-样本关联连线（最内圈）
for (i in seq_len(nrow(plot_data))) {
	circos.link(
		plot_data[i,2], c(accum_otu[plot_data[i,2]], accum_otu[plot_data[i,2]] - plot_data[i,4]),
		plot_data[i,1], c(accum_sample[plot_data[i,1]], accum_sample[plot_data[i,1]] - plot_data[i,3]),
		col = paste0(color_otu[plot_data[i,2]], '60'), border = NA )
	
	circos.rect(accum_otu[plot_data[i,2]], 0, accum_otu[plot_data[i,2]] - plot_data[i,4], 1, sector.index = plot_data[i,2], col = color_sample[plot_data[i,1]], border = NA)
	circos.rect(accum_sample[plot_data[i,1]], 0, accum_sample[plot_data[i,1]] - plot_data[i,3], 1, sector.index = plot_data[i,1], col = color_otu[plot_data[i,2]], border = NA)
	
	accum_otu[plot_data[i,2]] = accum_otu[plot_data[i,2]] - plot_data[i,4]
	accum_sample[plot_data[i,1]] = accum_sample[plot_data[i,1]] - plot_data[i,3]
}

##添加图例
otu_legend <- Legend(
		at = all_otu, labels = taxonomy$detail, labels_gp = gpar(fontsize = 8),    
		grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'), type = 'points', pch = NA, background = color_otu)

pushViewport(viewport(x = 0.85, y = 0.5))
grid.draw(otu_legend)
upViewport()
		
##清除 circlize 样式并关闭画板
circos.clear()
dev.off()



