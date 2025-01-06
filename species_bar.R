 

#color_its--------------------------------
cols1 <- c("#BEB8DC", "#FFBE7A")

palette <- c("#2878b5","#9ac9db","#ff8884", 
             "#4DBBD5E5","#00A087E5","#3C5488E5",
             "#F39B7FE5","#8491B4E5","#91D1C2E5",
             "#7E6148E5","#B09C85E5","#3B4992E5", 
             "#002EA6","#FFE76F","#01847F",
             "#008B45E5", "#631879E5", "#008280E5"  ,
             "#5F559BE5", "#00A087FF","#FF7802", "#91D1C2FF","#f8ac8c",
             "#3C5488FF","#000000","#BC3C29E5" ,"#0072B5E5" ,"#E18727E5",
             "#20854EE5" , "#7876B1E5" ,"#6F99ADE5" ,"#FFDC91E5", "#EE4C97E5")


pal <- c("#BEB8DC", "#FFBE7A",
         "#E64B35FF","#4DBBD5FF","#00A087FF","#F39B7FFF","#3C5488FF",
         "#E64B35FF","#4DBBD5FF","#00A087FF","#F39B7FFF","#3C5488FF",
         "#8491B4E5","#91D1C2E5","#7E6148E5",
         "#8491B4E5","#91D1C2E5","#7E6148E5",
         "#8491B4E5","#91D1C2E5","#7E6148E5",
         "#8491B4E5","#91D1C2E5","#7E6148E5",
         "#8491B4E5","#91D1C2E5","#7E6148E5",
         "#8491B4E5","#91D1C2E5","#7E6148E5",
         "#8491B4E5","#91D1C2E5","#7E6148E5",
         "#8491B4E5","#91D1C2E5","#7E6148E5",
         "#8491B4E5","#91D1C2E5","#7E6148E5",
         "#8491B4E5","#91D1C2E5","#7E6148E5")

#its------------------
taxa <- read_csv("../data/group_zgj.csv")
x1 <- read_tsv("../data/its_genus.tsv",col_names = TRUE) %>% t()
colnames(x1) <- x1[1,]
x1 <- x1[-1,] %>% as.data.frame()
x1 <- apply(x1, 2, as.numeric)
x1 <- x1 %>% as.data.frame() %>% mutate(Taxa=taxa$Sampleid) 
x1 <- x1 %>% remove_rownames()
#%>%   column_to_rownames("Taxa")
#apply(x1, 1, sum)

# x1 <- x1 %>% column_to_rownames("Taxa")
# rownames(x1) <- taxa$Sampleid

x1 <- x1 %>% pivot_longer(cols = !Taxa,names_to = "Samples",
             values_to = "number") 
colnames(x1)[1] <- "Sampleid"
#taxa <- taxa[-5]

meta_taxa <- taxa %>% 
  inner_join(.,x1,by="Sampleid")

title <- taxa$Sampleid
meta_taxa$Sampleid <- factor(meta_taxa$Sampleid , levels = title , ordered = T)
meta_taxa$Site <- factor(meta_taxa$Site ,
                                        levels = unique(meta_taxa$Site ) , ordered = T)
meta_taxa$Time <- factor(meta_taxa$Time ,
                         levels = unique(meta_taxa$Time ) , ordered = T)
meta_taxa$Repeat<- factor(meta_taxa$Repeat ,
                                      levels = c(1:10) , ordered = T)


its1 <- read_tsv("../data/its_genus.tsv",col_names = TRUE) %>% as.data.frame()
rownames(its1) <- its1[,1]
its1 <- its1[-1]
sum2 <- apply(its1, 1, sum) %>% as.data.frame()
orsum <- rownames(sum2)
orsum <- orsum[-14]
order<- c("Others",orsum)
meta_taxa$Samples <- factor(meta_taxa$Samples,levels = orsum,ordered = T)



p2 <- ggplot(meta_taxa,aes(Sampleid,number,fill=Samples))+
  geom_col(position="stack") +
  facet_nested(.~Site+Time+Repeat,drop=T,scale="free",
               space="free",switch="x",
               strip =strip_nested(
                 background_x = elem_list_rect(fill =pal),
                 by_layer_x = F
               ))+
  scale_fill_manual(values=palette,limits = order)+
  labs(x=NULL, y="Percent Genus Abundance")+
  scale_y_continuous(expand = c(0,0),labels=scales::percent)+
  theme(strip.background = element_rect(fill="white",color="black"),
        panel.spacing = unit(0,"lines"),
        strip.text.x = element_text(size=7,color="black"),
        axis.text.y=element_text(size=8,color="black"),
        axis.title.y = element_text(size=10,color="black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.key=element_blank(),   # 图例键为空
        #legend.position = "top",
        legend.text = element_text(color="black",size=9,face="italic"), # 定义图例文本
        legend.spacing.x=unit(0.1,'cm'), # 定义文本书平距离
        legend.key.width=unit(0.3,'cm'), # 定义图例水平大小
        legend.key.height=unit(0.3,'cm'), # 定义图例垂直大小
        legend.background=element_blank(),
        panel.grid.major=element_blank(), # 移除主网格线
        panel.grid.minor=element_blank(),
        legend.position = "bottom")+
  # guides(fill = guide_legend( ncol = 1, byrow = TRUE))+#图例设为i一列
  guides(fill=guide_legend(reverse = TRUE,
                           title = NULL,ncol = 7,
                           keywidth = unit(8,units = "points"),
                           keyheight = unit(8,units = "points")))
  
p2

ggsave(p2,device = "pdf",file="../result/species_bar/its_genus_bar_20_rev.pdf",width = 21, height = 16, units = "cm",onefile=F)   
ggsave(p2,device = "png",file="../result/species_bar/its_genus_bar_20_rev.png",
       dpi = 600,width = 24, height = 18, units = "cm") 


its <- read_tsv("../data/its_genus.tsv",col_names = TRUE) %>% t()
group <- read_csv("../data/group_zgj.csv")
colnames(its) <- its[1,]
its <-its[-1,] %>% as.data.frame()
its <- apply(its, 2, as.numeric)
its <- its %>% as.data.frame() %>% mutate(Sampleid=taxa$Sampleid) 

#colnames(its)[1] <- "Sampleid"
its1 <- its %>% left_join(group,by = "Sampleid")

its2 <- its1 %>% select(-Sampleid) %>% group_by(Site,Time) %>% 
  summarise(across(.fns = mean))

colnames(its2)[1] <- "Sampleid"
name <- c(paste0(its2$Sampleid,its2$Time)) %>% as.data.frame()
colnames(name)
its3 <- its2[,-c(1,2,24)] %>% mutate(Taxa = paste0(its2$Sampleid,its2$Time))
x1 <- its3 %>% pivot_longer(cols = !Taxa,names_to = "Samples",
                            values_to = "number") 
colnames(x1)[1] <- "Sampleid"
title1 <- c("CK0d", "CK3d","CK10d", "CK30d", "CK45d","Lpscw0d", "Lpscw3d","Lpscw10d", "Lpscw30d", "Lpscw45d")
x1$Sampleid <- factor(x1$Sampleid , levels = title1 , ordered = T)
x1$Samples <- factor(x1$Samples,levels = order,ordered = T)



q2 <- ggplot(x1,aes(Sampleid,number,fill=Samples))+
  geom_bar(stat = "identity",position = "stack",
           width = 0.8) +
  scale_fill_manual(values=palette,limits = order)+
  labs(x=NULL, y="Average Percent Genus Abundance")+
  scale_y_continuous(expand = c(0,0),labels=scales::percent)+
  theme(strip.background = element_rect(fill="white",color="black"),
        panel.spacing = unit(0,"lines"),
        strip.text.x = element_text(size=8,color="black"),
        axis.text.x=element_text(angle = 90,size=8,color="black",hjust = 1,vjust = 0.5),
        axis.text.y=element_text(size=8,color="black"),
        axis.title.y = element_text(size=8,color="black"),
        legend.key=element_blank(),   # 图例键为空
        #legend.position = "top",
        legend.text = element_text(color="black",size=9,face="italic"), # 定义图例文本
        legend.spacing.x=unit(0.1,'cm'), # 定义文本书平距离
        legend.key.width=unit(0.3,'cm'), # 定义图例水平大小
        legend.key.height=unit(0.3,'cm'), # 定义图例垂直大小
        legend.background=element_blank(),
        panel.grid.major=element_blank(), # 移除主网格线
        panel.grid.minor=element_blank())+
  # guides(fill = guide_legend( ncol = 1, byrow = TRUE))+#图例设为i一列
  guides(fill=guide_legend(reverse = TRUE,
                           title = NULL,ncol = 1,
                           keywidth = unit(8,units = "points"),
                           keyheight = unit(8,units = "points")))

q2

ggsave(q2,device = "pdf",file="../result/species_bar/its_genus_aev.pdf",width = 16, height = 16, units = "cm",onefile=F)   
ggsave(q2,device = "png",file="../result/species_bar/its_genus_aev.png",
       dpi = 300,width = 10, height = 10, units = "cm") 

ggsave(p2,device = "pdf",file="../result/species_bar/its_genus_bar_20.pdf",width = 21, height = 16, units = "cm",onefile=F)   
ggsave(p2,device = "png",file="../result/species_bar/its_genus_bar_20.png",
       dpi = 300,width = 14, height = 10, units = "cm") 


#color_16s--------------------------------
cols1 <- c("#BEB8DC", "#FFBE7A")

palette <- c("#2878b5","#9ac9db","#ff8884", 
             "#4DBBD5E5","#00A087E5","#3C5488E5",
             "#F39B7FE5","#8491B4E5","#91D1C2E5",
             "#7E6148E5","#B09C85E5","#3B4992E5", 
             "#008B45E5", "#631879E5", "#008280E5" ,
             "#5F559BE5", "#00A087FF",
             "#002EA6","#FFE76F","#01847F",
             "#E18727E5","#FF7802", "#91D1C2FF","#f8ac8c",
             "#3C5488FF","#000000","#BC3C29E5" ,"#0072B5E5" ,"#E18727E5",
             "#20854EE5" , "#7876B1E5" ,"#6F99ADE5" ,"#FFDC91E5", "#EE4C97E5")


pal <- c("#B883D4", "#F0988C",
         "#91D1C2FF","#4DBBD5E5","#00A087E5","#008B45E5","#20854EE5",
         "#91D1C2FF","#4DBBD5E5","#00A087E5","#008B45E5","#20854EE5",
         "#FFDC91E5","#f8ac8c","#E18727E5",
         "#FFDC91E5","#f8ac8c","#E18727E5",
         "#FFDC91E5","#f8ac8c","#E18727E5",
         "#FFDC91E5","#f8ac8c","#E18727E5",
         "#FFDC91E5","#f8ac8c","#E18727E5",
         "#FFDC91E5","#f8ac8c","#E18727E5",
         "#FFDC91E5","#f8ac8c","#E18727E5",
         "#FFDC91E5","#f8ac8c","#E18727E5",
         "#FFDC91E5","#f8ac8c","#E18727E5",
         "#FFDC91E5","#f8ac8c","#E18727E5")




#16s-----------------------------------------------
taxa <- read_csv("../data/group_zgj.csv")
x1 <- read_tsv("../data/16s_genus.tsv",col_names = TRUE) %>% t()
colnames(x1) <- x1[1,]
x1 <- x1[-1,] %>% as.data.frame()
x1 <- apply(x1, 2, as.numeric)
x1 <- x1 %>% as.data.frame() %>% mutate(Taxa=taxa$Sampleid) 
x1 <- x1 %>% remove_rownames()
#%>%   column_to_rownames("Taxa")
#apply(x1, 1, sum)

# x1 <- x1 %>% column_to_rownames("Taxa")
# rownames(x1) <- taxa$Sampleid

x1 <- x1 %>% pivot_longer(cols = !Taxa,names_to = "Samples",
                          values_to = "number") 
colnames(x1)[1] <- "Sampleid"
#taxa <- taxa[-5]

meta_taxa <- taxa %>% 
  inner_join(.,x1,by="Sampleid")

title <- taxa$Sampleid
meta_taxa$Sampleid <- factor(meta_taxa$Sampleid , levels = title , ordered = T)
meta_taxa$Site <- factor(meta_taxa$Site ,
                         levels = unique(meta_taxa$Site ) , ordered = T)
meta_taxa$Time <- factor(meta_taxa$Time ,
                         levels = unique(meta_taxa$Time ) , ordered = T)
meta_taxa$Repeat<- factor(meta_taxa$Repeat ,
                          levels = c(1:10) , ordered = T)


s161 <- read_tsv("../data/16s_genus.tsv",col_names = TRUE) %>% as.data.frame()
rownames(s161) <- s161[,1]
s161 <- s161[-1]
sum2 <- apply(s161, 1, sum) %>% as.data.frame()
orsum <- rownames(sum2)
orsum <- orsum[-17]
order<- c("Others",orsum)
meta_taxa$Samples <- factor(meta_taxa$Samples,levels = order,ordered = T)



p4 <- ggplot(meta_taxa,aes(Sampleid,number,fill=Samples))+
  geom_col(position="stack") +
  facet_nested(.~Site+Time+Repeat,drop=T,scale="free",
               space="free",switch="x",
               strip =strip_nested(
                 background_x = elem_list_rect(fill =pal),
                 by_layer_x = F
               ))+
  scale_fill_manual(values=palette,limits = order)+
  labs(x=NULL, y="Percent Genus Abundance")+
  scale_y_continuous(expand = c(0,0),labels=scales::percent)+
  theme(strip.background = element_rect(fill="white",color="black"),
        panel.spacing = unit(0,"lines"),
        strip.text.x = element_text(size=7,color="black"),
        axis.text.y=element_text(size=8,color="black"),
        axis.title.y = element_text(size=10,color="black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.key=element_blank(),   # 图例键为空
        #legend.position = "top",
        legend.text = element_text(color="black",size=9,face="italic"), # 定义图例文本
        legend.spacing.x=unit(0.1,'cm'), # 定义文本书平距离
        legend.key.width=unit(0.3,'cm'), # 定义图例水平大小
        legend.key.height=unit(0.3,'cm'), # 定义图例垂直大小
        legend.background=element_blank(),
        panel.grid.major=element_blank(), # 移除主网格线
        panel.grid.minor=element_blank(),
        legend.position = "bottom")+
  # guides(fill = guide_legend( ncol = 1, byrow = TRUE))+#图例设为i一列
  guides(fill=guide_legend(reverse = TRUE,
                           title = NULL,ncol = 7,
                           keywidth = unit(8,units = "points"),
                           keyheight = unit(8,units = "points")))

p4

ggsave(p4,device = "pdf",file="../result/species_bar/16s_genus_bar_20_rev.pdf",width = 21, height = 16, units = "cm",onefile=F)   
ggsave(p4,device = "png",file="../result/species_bar/16s_genus_bar_20_rev.png",
       dpi = 600,width = 26, height = 19, units = "cm") 



its <- read_tsv("../data/16s_genus.tsv",col_names = TRUE) %>% t()
group <- read_csv("../data/group_zgj.csv")
colnames(its) <- its[1,]
its <-its[-1,] %>% as.data.frame()
its <- apply(its, 2, as.numeric)
its <- its %>% as.data.frame() %>% mutate(Sampleid=taxa$Sampleid) 

#colnames(its)[1] <- "Sampleid"
its1 <- its %>% left_join(group,by = "Sampleid")

its2 <- its1 %>% select(-Sampleid) %>% group_by(Site,Time) %>% 
  summarise(across(.fns = mean))

colnames(its2)[1] <- "Sampleid"
name <- c(paste0(its2$Sampleid,its2$Time)) %>% as.data.frame()
colnames(name)
its3 <- its2[,-c(1,2,24)] %>% mutate(Taxa = paste0(its2$Sampleid,its2$Time))
x1 <- its3 %>% pivot_longer(cols = !Taxa,names_to = "Samples",
                            values_to = "number") 
colnames(x1)[1] <- "Sampleid"

title1 <- c("CK0d", "CK3d","CK10d", "CK30d", "CK45d","Lpscw0d", "Lpscw3d","Lpscw10d", "Lpscw30d", "Lpscw45d")
x1$Sampleid <- factor(x1$Sampleid , levels = title1 , ordered = T)
x1$Samples <- factor(x1$Samples,levels = order,ordered = T)



q4 <- ggplot(x1,aes(Sampleid,number,fill=Samples))+
  geom_bar(stat = "identity",position = "stack",
           width = 0.8) +
  scale_fill_manual(values=palette,limits = order)+
  labs(x=NULL, y="Average Percent Genus Abundance")+
  scale_y_continuous(expand = c(0,0),labels=scales::percent)+
  theme(strip.background = element_rect(fill="white",color="black"),
        panel.spacing = unit(0,"lines"),
        strip.text.x = element_text(size=8,color="black"),
        axis.text.x=element_text(angle = 90,size=8,color="black",hjust = 1,vjust = 0.5),
        axis.text.y=element_text(size=8,color="black"),
        axis.title.y = element_text(size=8,color="black"),
        legend.key=element_blank(),   # 图例键为空
        #legend.position = "top",
        legend.text = element_text(color="black",size=9,face="italic"), # 定义图例文本
        legend.spacing.x=unit(0.1,'cm'), # 定义文本书平距离
        legend.key.width=unit(0.3,'cm'), # 定义图例水平大小
        legend.key.height=unit(0.3,'cm'), # 定义图例垂直大小
        legend.background=element_blank(),
        panel.grid.major=element_blank(), # 移除主网格线
        panel.grid.minor=element_blank())+
  # guides(fill = guide_legend( ncol = 1, byrow = TRUE))+#图例设为i一列
  guides(fill=guide_legend(reverse = TRUE,
                           title = NULL,ncol = 1,
                           keywidth = unit(8,units = "points"),
                           keyheight = unit(8,units = "points")))

q4

ggsave(q4,device = "pdf",file="../result/species_bar/16s_genus_aev.pdf",width = 16, height = 16, units = "cm",onefile=F)   
ggsave(q4,device = "png",file="../result/species_bar/16s_genus_aev.png",
       dpi = 300,width = 10, height = 10, units = "cm") 

ggsave(p4,device = "pdf",file="../result/species_bar/16s_genus_bar_20.pdf",width = 21, height = 16, units = "cm",onefile=F)   
ggsave(p4,device = "png",file="../result/species_bar/16s_genus_bar_20.png",
       dpi = 300,width = 14, height = 10, units = "cm") 
