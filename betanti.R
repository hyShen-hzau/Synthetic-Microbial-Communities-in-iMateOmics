library(tidyverse)
library(vegan)
group <- read_tsv("metadata.txt", comment = "#")
#group_with_starter <- read_tsv("metadata_with_starter.txt", comment = "#")

#16S

bac_bnti <- read.delim("16s/betanti/bacteria_betanti_raxml.tsv",row.names = 1)
# bac_bnti <- read.delim("16s/betanti/bac_betanti.tsv",row.names = 1)

diag(bac_bnti) <- NA
bac_bnti[upper.tri(bac_bnti)] <- NA
bac_bnti_dist <- bac_bnti

bac_bnti <- bac_bnti %>% rownames_to_column("sample1") %>% 
  pivot_longer(cols = -1,names_to="sample2",values_to="betanti") %>% filter(!is.na(betanti))

bac_bnti <- left_join(bac_bnti,group,by=c("sample1"="sampleid")) %>% left_join(group,by=c("sample2"="sampleid"),suffix=c("1","2"))


# whole
bac_bnti %>% 
  filter(SampleTime1==SampleTime2,StarterType1==StarterType2) %>% 
  ggplot(aes(SampleTime1,betanti))+
  geom_boxplot()+geom_jitter()+
  geom_hline(yintercept = c(-2,2),linetype="dashed")+
  scale_x_discrete(limits=unique(group$SampleTime))+
  theme_bw()+theme(panel.grid = element_blank())
  
# pure

bac_bnti %>% 
  filter(StarterType1=="PureStarter"&StarterType2=="PureStarter") %>% 
  filter(SampleTime1==SampleTime2) %>% 
  ggplot(aes(SampleTime1,betanti))+
  geom_boxplot()+geom_jitter()+
  geom_hline(yintercept = c(-2,2),linetype="dashed")+
  scale_x_discrete(limits=unique(group$SampleTime))+
  theme_bw()+theme(panel.grid = element_blank())

#traditional  

bac_bnti %>% 
  filter(StarterType1=="TraditionalStarter"&StarterType2=="TraditionalStarter") %>% 
  filter(SampleTime1==SampleTime2) %>% 
  ggplot(aes(SampleTime1,betanti))+
  geom_boxplot()+geom_jitter()+
  geom_hline(yintercept = c(-2,2),linetype="dashed")+
  scale_x_discrete(limits=unique(group$SampleTime))+
  theme_bw()+theme(panel.grid = element_blank())


# lm

bac_bnti <- bac_bnti %>% filter(SampleTime1==SampleTime2,StarterType1==StarterType2)
  
bac_bnti$SampleTime1 <- str_replace_all(bac_bnti$SampleTime1,
              c("SAC0h"="1","SAC6h"="2","SAC12h"="3",
                "SAC18h"="4","SAC24h"="5",
                "FER0d"="6","FER3d"="7","FER6d"="8",
                "FER9d"="9","FER12d"="10","FER15d"="11"))
bac_bnti$SampleTime1 <- as.numeric(bac_bnti$SampleTime1)

library(segmented)
bac_lmmodel <- lm(betanti~SampleTime1,data = bac_bnti)

summary(lm(betanti~SampleTime1,data = bac_bnti))
ggplot(bac_bnti,aes(SampleTime1,betanti))+geom_smooth(method = "lm")+geom_point()

bac_segmodel <- segmented(bac_lmmodel,seg.Z = ~SampleTime1,npsi=2)
plot(bac_segmodel)
summary(bac_segmodel)
# result
# 4.000 8.000
# phase 1 -1.0197x+0.3600


bac_bnti1 <- bac_bnti %>% filter(SampleTime1 %in% c(1:4))
ggplot(bac_bnti1,aes(SampleTime1,betanti))+geom_smooth(method = "lm")
summary(lm(betanti~SampleTime1,bac_bnti1))

bac_bnti2 <- bac_bnti %>% filter(SampleTime1 %in% c(4:7))
ggplot(bac_bnti2,aes(SampleTime1,betanti))+geom_smooth(method = "lm")
summary(lm(betanti~SampleTime1,bac_bnti2))

bac_bnti3 <- bac_bnti %>% filter(SampleTime1 %in% c(7:11))
ggplot(bac_bnti3,aes(SampleTime1,betanti))+geom_smooth(method = "lm")
summary(lm(betanti~SampleTime1,bac_bnti3))

bac_bnti %>%
  ggplot(aes(SampleTime1,betanti,group=SampleTime1))+
  geom_boxplot(alpha=0.5,outlier.colour = NA)+geom_jitter(alpha=0.5)+
  geom_hline(yintercept = c(-2,2),linetype="dashed")+
  theme_bw()+theme(panel.grid = element_blank(),axis.title.x = element_blank())+
  scale_x_discrete(limits=unique(group$SampleTime))+
  scale_y_continuous(breaks = c(-4,-2,0,2,4,8))+
  labs(y="βNTI")+
  geom_smooth(aes(SampleTime1,betanti,group=1),color="#E41A1C",method = "lm",data = bac_bnti1,se=F)+
  geom_smooth(aes(SampleTime1,betanti,group=1),color="#377EB8",method = "lm",data = bac_bnti2,se=F)+
  geom_smooth(aes(SampleTime1,betanti,group=1),color="#4DAF4A",method = "lm",data = bac_bnti3,se=F)

# bac_bnti4 <- bac_bnti %>% filter(SampleTime1 %in% c(1:7))
# bac_bnti5 <- bac_bnti %>% filter(SampleTime1 %in% c(8:11))
# 
# loess(betanti~SampleTime1,bac_bnti)
# 
# bac_bnti %>%
#   ggplot(aes(SampleTime1,betanti,group=SampleTime1))+
#   geom_boxplot(alpha=0.5,outlier.color = NA)+geom_jitter(alpha=0.5)+
#   geom_smooth(aes(group=1),se=F,data = bac_bnti4)+
#   geom_smooth(aes(group=1),se=F,data = bac_bnti5)+
#   geom_hline(yintercept = c(-2,2),linetype="dashed")+
#     theme_bw()+theme(panel.grid = element_blank(),axis.title.x = element_blank())+
#     scale_x_discrete(limits=unique(group$SampleTime))+
#     scale_y_continuous(breaks = c(-4,-2,0,2,4,8))+
#     labs(y="βNTI")

bac_bnti %>%
  ggplot(aes(SampleTime1,betanti,group=SampleTime1))+
  geom_boxplot(outlier.color = NA)+geom_jitter(aes(color=StarterType1))+
  geom_smooth(aes(group=1),se=F,data = bac_bnti)+
  geom_hline(yintercept = c(-2,2),linetype="dashed")+
  theme_bw()+theme(panel.grid = element_blank(),axis.title.x = element_blank(),
                   axis.text.x = element_text(angle = 90,vjust = 0.5))+
  scale_x_discrete(limits=unique(group$SampleTime))+
  scale_y_continuous(breaks = c(-4,-2,0,2,4,8))+
  labs(y="βNTI")+
  scale_color_manual(values=c("PureStarter"="#00BFC4","TraditionalStarter"="#F8766D"))+
  labs(color="Batch")

#itS

its_bnti <- read.delim("its/betanti/fungi_betaNTI.tsv",row.names = 1)

diag(its_bnti) <- NA
its_bnti[upper.tri(its_bnti)] <- NA

its_bnti <- its_bnti %>% rownames_to_column("sample1")
its_bnti <- its_bnti %>% pivot_longer(cols = -1,names_to="sample2",values_to="betanti") %>% filter(!is.na(betanti))

its_bnti <- left_join(its_bnti,group,by=c("sample1"="sampleid")) %>% left_join(group,by=c("sample2"="sampleid"),suffix=c("1","2"))

# whole
its_bnti %>% 
  filter(SampleTime1==SampleTime2) %>% 
  ggplot(aes(SampleTime1,betanti))+
  geom_boxplot()+geom_jitter()+
  geom_hline(yintercept = c(-2,2),linetype="dashed")+
  scale_x_discrete(limits=unique(group$SampleTime))+
  theme_bw()+theme(panel.grid = element_blank())

# pure

its_bnti %>% 
  filter(StarterType1=="PureStarter"&StarterType2=="PureStarter") %>% 
  filter(SampleTime1==SampleTime2) %>% 
  ggplot(aes(SampleTime1,betanti))+
  geom_boxplot()+geom_jitter()+
  geom_hline(yintercept = c(-2,2),linetype="dashed")+
  scale_x_discrete(limits=unique(group$SampleTime))+
  theme_bw()+theme(panel.grid = element_blank())

#traditional  

  its_bnti %>% 
  filter(StarterType1=="TraditionalStarter"&StarterType2=="TraditionalStarter") %>% 
  filter(SampleTime1==SampleTime2) %>% 
  ggplot(aes(SampleTime1,betanti))+
  geom_boxplot()+geom_jitter()+
  geom_hline(yintercept = c(-2,2),linetype="dashed")+
  scale_x_discrete(limits=unique(group$SampleTime))+
  theme_bw()+theme(panel.grid = element_blank())

its_bnti <- its_bnti %>% filter(SampleTime1==SampleTime2)

its_bnti$SampleTime1 <- str_replace_all(its_bnti$SampleTime1,
                                        c("SAC0h"="1","SAC6h"="2","SAC12h"="3",
                                          "SAC18h"="4","SAC24h"="5",
                                          "FER0d"="6","FER3d"="7","FER6d"="8",
                                          "FER9d"="9","FER12d"="10","FER15d"="11"))
its_bnti$SampleTime1 <- as.numeric(its_bnti$SampleTime1)


its_lmmodel <- lm(betanti~SampleTime1,data = its_bnti)
summary(lm(betanti~SampleTime1,data = its_bnti))
ggplot(its_bnti,aes(SampleTime1,betanti))+geom_smooth(method = "lm")+geom_point()

its_bnti %>% 
  ggplot(aes(SampleTime1,betanti,group=SampleTime1))+
  geom_boxplot(outlier.colour = NA)+geom_jitter(aes(color=StarterType1))+
  geom_hline(yintercept = c(-2,2),linetype="dashed")+
  scale_x_discrete(limits=unique(group$SampleTime))+
  labs(y="βNTI")+
  theme_bw()+theme(panel.grid = element_blank(),axis.title.x = element_blank(),
                   axis.text.x = element_text(angle = 90,vjust = 0.5))+
  geom_smooth(aes(group=1),color="#984EA3",method = "lm",se=F)+
  scale_color_manual(values=c("PureStarter"="#00BFC4","TraditionalStarter"="#F8766D"))+
  labs(color="Batch")

# env 

env <- read_tsv("env.txt")
env <- left_join(env,group)
env <- env %>% column_to_rownames("sampleid")
env <- env[rownames(bac_bnti_dist),] %>% rownames_to_column("sampleid")

env1 <- env %>% filter(SampleTime %in% c("SAC0h","SAC6h","SAC12h","SAC18h"))
env3 <- env %>% filter(SampleTime %in% c("FER3d","FER6d","FER9d","FER12d","FER15d"))
# env3 <- env %>% filter(SampleTime %in% c("FER0d","FER3d"))

# phase 1

bac_bnti1_dis <- bac_bnti_dist[c(1:12,34:45),c(1:12,34:45)]
bac_bnti1_dis <- bac_bnti1_dis %>% rownames_to_column("sample1")
bac_bnti1_dis <- bac_bnti1_dis %>% pivot_longer(cols = -1,names_to="sample2",values_to="betanti") %>% filter(!is.na(betanti))
bac_bnti1_dis <- left_join(bac_bnti1_dis,group,by=c("sample1"="sampleid")) %>% left_join(group,by=c("sample2"="sampleid"),suffix=c("1","2"))
bac_bnti1_dis <- filter(bac_bnti1_dis,StarterType1==StarterType2,SampleTime1==SampleTime2) %>% select(SampleTime=SampleTime1,betanti)

# Moisture 

Moisture1 <- env1 %>% select(sampleid,Moisture) %>% left_join(group) %>% select(SampleTime,Moisture)

bnti_Moisture1 <- left_join(bac_bnti1_dis,Moisture1)
summary(lm(betanti~log(Moisture),data = bnti_Moisture1))
ggplot(bnti_Moisture1,aes(log(Moisture),betanti))+geom_point()+geom_smooth(method = "lm")

# ReducingSugar

ReducingSugar1 <- env1 %>% select(sampleid,ReducingSugar) %>% left_join(group) %>% select(SampleTime,ReducingSugar)

bnti_ReducingSugar1 <- left_join(bac_bnti1_dis,ReducingSugar1)
summary(lm(betanti~log(ReducingSugar),data = bnti_ReducingSugar1))
ggplot(bnti_ReducingSugar1,aes(log(ReducingSugar),betanti))+geom_jitter()+geom_smooth(method = "lm",se=F)+
  theme_bw(base_size = 8)+
  theme(panel.grid = element_blank())+
  geom_hline(yintercept = c(-2,2),linetype="dashed")+
  scale_y_continuous(breaks=c(-5,-4,-3,-2,-1,0,1,2))

# pdf("../paper/prepare/phase1_reducingSugar.pdf",width = 55.866/25.4,height = 46.019/25.4,family = "ArialMT",version = "1.5",useDingbats = F)
# ggplot(bnti_ReducingSugar1,aes(log(ReducingSugar),betanti))+geom_jitter()+geom_smooth(method = "lm",se=F)+
#   theme_bw(base_size = 8)+
#   theme(panel.grid = element_blank())+
#   geom_hline(yintercept = c(-2,2),linetype="dashed")+
#   scale_y_continuous(breaks=c(-5,-4,-3,-2,-1,0,1,2))
# dev.off()
# Acidity

Acidity1 <- env1 %>% select(sampleid,Acidity) %>% left_join(group) %>% select(SampleTime,Acidity)

bnti_Acidity1 <- left_join(bac_bnti1_dis,Acidity1)
summary(lm(betanti~log(Acidity),data = bnti_Acidity1))
ggplot(bnti_Acidity1,aes(log(Acidity),betanti))+geom_point()+geom_smooth(method = "lm")

# pH

pH1 <- env1 %>% select(sampleid,pH) %>% left_join(group) %>% select(SampleTime,pH)

bnti_pH1 <- left_join(bac_bnti1_dis,pH1)
summary(lm(betanti~pH,data = bnti_pH1))
ggplot(bnti_pH1,aes(pH,betanti))+geom_jitter()+geom_smooth(method = "lm",se=F)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_hline(yintercept = c(-2,2),linetype="dashed")+
  scale_y_continuous(breaks=c(-5,-4,-3,-2,-1,0,1,2))

# Temperature

Temperature1 <- env1 %>% select(sampleid,Temperature) %>% left_join(group) %>% select(SampleTime,Temperature)

bnti_Temperature1 <- left_join(bac_bnti1_dis,Temperature1)
summary(lm(betanti~Temperature,data = bnti_Temperature1))
ggplot(bnti_Temperature1,aes(Temperature,betanti))+geom_jitter(size=1)+geom_smooth(method = "lm",se=F)+
  theme_bw(base_size = 8)+
  theme(panel.grid = element_blank())+
  geom_hline(yintercept = c(-2,2),linetype="dashed")+
  scale_y_continuous(breaks=c(-5,-4,-3,-2,-1,0,1,2))+
  labs(x="Temperature( )",y="βNTI")

# Ethanol

Ethanol1 <- env1 %>% select(sampleid,Ethanol) %>% left_join(group) %>% select(SampleTime,Ethanol)

bnti_Ethanol1 <- left_join(bac_bnti1_dis,Ethanol1)
summary(lm(betanti~log(Ethanol),data = bnti_Ethanol1))
ggplot(bnti_Ethanol1,aes(log(Ethanol),betanti))+geom_point()+geom_smooth(method = "lm")

# phase 3

bac_bnti3_dis <- bac_bnti_dist[c(19:33,52:66),c(19:33,52:66)]
bac_bnti3_dis <- bac_bnti3_dis %>% rownames_to_column("sample1")
bac_bnti3_dis <- bac_bnti3_dis %>% pivot_longer(cols = -1,names_to="sample2",values_to="betanti") %>% filter(!is.na(betanti))
bac_bnti3_dis <- left_join(bac_bnti3_dis,group,by=c("sample1"="sampleid")) %>% left_join(group,by=c("sample1"="sampleid"),suffix=c("1","2"))
bac_bnti3_dis <- filter(bac_bnti3_dis,StarterType1==StarterType2,SampleTime1==SampleTime2) %>% select(SampleTime=SampleTime1,betanti)

# Moisture

Moisture3 <- env3 %>% select(sampleid,Moisture) %>% left_join(group) %>% select(SampleTime,Moisture)

bnti_Moisture3 <- left_join(bac_bnti3_dis,Moisture3)
summary(lm(betanti~log(Moisture),data = bnti_Moisture3))
ggplot(bnti_Moisture3,aes(log(Moisture),betanti))+geom_point()+geom_smooth(method = "lm")

# ReducingSugar

ReducingSugar3 <- env3 %>% select(sampleid,ReducingSugar) %>% left_join(group) %>% select(SampleTime,ReducingSugar)

bnti_ReducingSugar3 <- left_join(bac_bnti3_dis,ReducingSugar3)
summary(lm(betanti~log(ReducingSugar),data = bnti_ReducingSugar3))
ggplot(bnti_ReducingSugar3,aes(log(ReducingSugar),betanti))+geom_point()+geom_smooth(method = "lm")

# Acidity

Acidity3 <- env3 %>% select(sampleid,Acidity) %>% left_join(group) %>% select(SampleTime,Acidity)

bnti_Acidity3 <- left_join(bac_bnti3_dis,Acidity3)
summary(lm(betanti~log(Acidity),data = bnti_Acidity3))
ggplot(bnti_Acidity3,aes(log(Acidity),betanti))+geom_jitter()+geom_smooth(method = "lm")+
  theme_bw()+theme(panel.grid = element_blank())+
  geom_hline(yintercept = c(-2,2),linetype="dashed")+
  scale_y_continuous(breaks=c(-5,-2,0,2,5,10))

# pH

pH3 <- env3 %>% select(sampleid,pH) %>% left_join(group) %>% select(SampleTime,pH)

bnti_pH3 <- left_join(bac_bnti3_dis,pH3)
summary(lm(betanti~pH,data = bnti_pH3))
ggplot(bnti_pH3,aes(pH,betanti))+geom_point()+geom_smooth(method = "lm")

# Temperature

Temperature3 <- env3 %>% select(sampleid,Temperature) %>% left_join(group) %>% select(SampleTime,Temperature)

bnti_Temperature3 <- left_join(bac_bnti3_dis,Temperature3)
summary(lm(betanti~Temperature,data = bnti_Temperature3))
ggplot(bnti_Temperature3,aes(Temperature,betanti))+geom_jitter(size=1)+geom_smooth(method = "lm",se=F)+
  theme_bw(base_size = 8)+
  theme(panel.grid = element_blank())+
  geom_hline(yintercept = c(-2,2),linetype="dashed")+
  scale_y_continuous(breaks=c(-5,-2,0,2,5,10))+
  labs(x="Temperature( )",y="βNTI")

# Ethanol

Ethanol3 <- env3 %>% select(sampleid,Ethanol) %>% left_join(group) %>% select(SampleTime,Ethanol)

bnti_Ethanol3 <- left_join(bac_bnti3_dis,Ethanol3)
summary(lm(betanti~log(Ethanol),data = bnti_Ethanol3))
ggplot(bnti_Ethanol3,aes(log(Ethanol),betanti))+geom_point()+geom_smooth(method = "lm")



# env1 <- env %>% filter(SampleTime %in% c("SAC0h","SAC6h","SAC12h","SAC18h"))
# env3 <- env %>% filter(SampleTime %in% c("FER6d","FER9d","FER12d","FER15d"))
# 
# 
# # phase 1
# 
# bac_bnti1_dis <- bac_bnti_dist[c(1:12,34:45),c(1:12,34:45)]
# bac_bnti1_dis <- bac_bnti1_dis %>% rownames_to_column("sample1")
# bac_bnti1_dis <- bac_bnti1_dis %>% pivot_longer(cols = -1,names_to="sample2",values_to="betanti") %>% filter(!is.na(betanti))
# bac_bnti1_dis <- left_join(bac_bnti1_dis,group,by=c("sample1"="sampleid")) %>% left_join(group,by=c("sample2"="sampleid"),suffix=c("1","2"))
# bac_bnti1_dis <- filter(bac_bnti1_dis,StarterType1==StarterType2)
# bac_bnti1_dist <- bac_bnti1_dis[1:3] %>%
#   pivot_wider(names_from = "sample2",values_from="betanti") %>%
#   column_to_rownames("sample1") %>% as.dist
# 
# # Moisture
# Moisture1 <- env1 %>% select(sampleid,Moisture) %>% column_to_rownames("sampleid")
# Moisture1_dis <- vegdist(Moisture1,method = "euclidean",diag = T)
# 
# Moisture1_dis <- Moisture1_dis %>% as.matrix %>% data.frame(check.names = F,stringsAsFactors = F)
# diag(Moisture1_dis) <- NA ;Moisture1_dis[upper.tri(Moisture1_dis)] <- NA
# Moisture1_dis <- Moisture1_dis %>% rownames_to_column("sample1")
# Moisture1_dis <- Moisture1_dis %>% pivot_longer(cols = -1,names_to="sample2",values_to="Moisture") %>% filter(!is.na(Moisture))
# Moisture1_dis <- left_join(Moisture1_dis,group,by=c("sample1"="sampleid")) %>% left_join(group,by=c("sample2"="sampleid"),suffix=c("1","2"))
# Moisture1_dis <- filter(Moisture1_dis,StarterType1==StarterType2)
# Moisture1_dist <- Moisture1_dis[1:3] %>%
#   pivot_wider(names_from = "sample2",values_from="Moisture") %>%
#   column_to_rownames("sample1") %>% as.dist
# 
# bnti_Moisture1 <- left_join(bac_bnti1_dis[1:3],Moisture1_dis)
# summary(lm(betanti~Moisture,data = bnti_Moisture1))
# ggplot(bnti_Moisture1,aes(Moisture,betanti))+geom_point()+geom_smooth(method = "lm")
# 
# mantel(Moisture1_dist,bac_bnti1_dist,na.rm=T)
# 
# # Reducing Sugar
# ReducingSugar1 <- env1 %>% select(sampleid,ReducingSugar) %>% column_to_rownames("sampleid")
# ReducingSugar1_dis <- vegdist(ReducingSugar1,method = "euclidean",diag = T)
# 
# ReducingSugar1_dis <- ReducingSugar1_dis %>% as.matrix %>% data.frame(check.names = F,stringsAsFactors = F)
# diag(ReducingSugar1_dis) <- NA ;ReducingSugar1_dis[upper.tri(ReducingSugar1_dis)] <- NA
# ReducingSugar1_dis <- ReducingSugar1_dis %>% rownames_to_column("sample1")
# ReducingSugar1_dis <- ReducingSugar1_dis %>% pivot_longer(cols = -1,names_to="sample2",values_to="ReducingSugar") %>% filter(!is.na(ReducingSugar))
# ReducingSugar1_dis <- left_join(ReducingSugar1_dis,group,by=c("sample1"="sampleid")) %>% left_join(group,by=c("sample2"="sampleid"),suffix=c("1","2"))
# ReducingSugar1_dis <- filter(ReducingSugar1_dis,StarterType1==StarterType2)
# ReducingSugar1_dist <- ReducingSugar1_dis[1:3] %>%
#   pivot_wider(names_from = "sample2",values_from="ReducingSugar") %>%
#   column_to_rownames("sample1") %>% as.dist
# 
# 
# bnti_ReducingSugar1 <- left_join(bac_bnti1_dis[1:3],ReducingSugar1_dis)
# summary(lm(betanti~ReducingSugar,data = bnti_ReducingSugar1))
# ggplot(bnti_ReducingSugar1,aes(ReducingSugar,betanti))+geom_point()+geom_smooth(method = "lm")
# 
# mantel(ReducingSugar1_dist,bac_bnti1_dist,na.rm=T)
# 
# # Acidity
# 
# Acidity1 <- env1 %>% select(sampleid,Acidity) %>% column_to_rownames("sampleid")
# Acidity1_dis <- vegdist(Acidity1,method = "euclidean",diag = T)
# Acidity1_dis <- Acidity1_dis %>% as.matrix %>% data.frame(check.names = F,stringsAsFactors = F)
# diag(Acidity1_dis) <- NA ;Acidity1_dis[upper.tri(Acidity1_dis)] <- NA
# Acidity1_dis <- Acidity1_dis %>% rownames_to_column("sample1")
# Acidity1_dis <- Acidity1_dis %>% pivot_longer(cols = -1,names_to="sample2",values_to="Acidity") %>% filter(!is.na(Acidity))
# Acidity1_dis <- left_join(Acidity1_dis,group,by=c("sample1"="sampleid")) %>% left_join(group,by=c("sample2"="sampleid"),suffix=c("1","2"))
# Acidity1_dis <- filter(Acidity1_dis,StarterType1==StarterType2)
# Acidity1_dist <- Acidity1_dis[1:3] %>%
#   pivot_wider(names_from = "sample2",values_from="Acidity") %>%
#   column_to_rownames("sample1") %>% as.dist
# 
# bnti_Acidity1 <- left_join(bac_bnti1_dis[1:3],Acidity1_dis)
# summary(lm(betanti~Acidity,data = bnti_Acidity1))
# ggplot(bnti_Acidity1,aes(Acidity,betanti))+geom_point()+geom_smooth(method = "lm")
# 
# mantel(Acidity1_dist,bac_bnti1_dist,na.rm=T)
# 
# # Temperature
# Temperature1 <- env1 %>% select(sampleid,Temperature) %>% column_to_rownames("sampleid")
# Temperature1_dis <- vegdist(Temperature1,method = "euclidean",diag = T)
# Temperature1_dis <- Temperature1_dis %>% as.matrix %>% data.frame(check.names = F,stringsAsFactors = F)
# diag(Temperature1_dis) <- NA ;Temperature1_dis[upper.tri(Temperature1_dis)] <- NA
# Temperature1_dis <- Temperature1_dis %>% rownames_to_column("sample1")
# Temperature1_dis <- Temperature1_dis %>% pivot_longer(cols = -1,names_to="sample2",values_to="Temperature") %>% filter(!is.na(Temperature))
# Temperature1_dis <- left_join(Temperature1_dis,group,by=c("sample1"="sampleid")) %>% left_join(group,by=c("sample2"="sampleid"),suffix=c("1","2"))
# Temperature1_dis <- filter(Temperature1_dis,StarterType1==StarterType2)
# Temperature1_dist <- Temperature1_dis[1:3] %>%
#   pivot_wider(names_from = "sample2",values_from="Temperature") %>%
#   column_to_rownames("sample1") %>% as.dist
# 
# bnti_Temperature1 <- left_join(bac_bnti1_dis[1:3],Temperature1_dis)
# summary(lm(betanti~Temperature,data = bnti_Temperature1))
# ggplot(bnti_Temperature1,aes(Temperature,betanti))+geom_point()+geom_smooth(method = "lm")
# 
# mantel(Temperature1_dist,bac_bnti1_dist,na.rm=T)
# 
# # pH
# 
# pH1 <- env1 %>% select(sampleid,pH) %>% column_to_rownames("sampleid")
# pH1_dis <- vegdist(pH1,method = "euclidean",diag = T)
# pH1_dis <- pH1_dis %>% as.matrix %>% data.frame(check.names = F,stringsAsFactors = F)
# diag(pH1_dis) <- NA ;pH1_dis[upper.tri(pH1_dis)] <- NA
# pH1_dis <- pH1_dis %>% rownames_to_column("sample1")
# pH1_dis <- pH1_dis %>% pivot_longer(cols = -1,names_to="sample2",values_to="pH") %>% filter(!is.na(pH))
# pH1_dis <- left_join(pH1_dis,group,by=c("sample1"="sampleid")) %>% left_join(group,by=c("sample2"="sampleid"),suffix=c("1","2"))
# pH1_dis <- filter(pH1_dis,StarterType1==StarterType2)
# pH1_dist <- pH1_dis[1:3] %>%
#   pivot_wider(names_from = "sample2",values_from="pH") %>%
#   column_to_rownames("sample1") %>% as.dist
# 
# bnti_pH1 <- left_join(bac_bnti1_dis[1:3],pH1_dis)
# summary(lm(betanti~pH,data = bnti_pH1))
# ggplot(bnti_pH1,aes(pH,betanti))+geom_point()+geom_smooth(method = "lm")
# 
# mantel(pH1_dist,bac_bnti1_dist,na.rm = T)
# 
# # Ethanol
# Ethanol1 <- env1 %>% select(sampleid,Ethanol) %>% column_to_rownames("sampleid")
# Ethanol1_dis <- vegdist(Ethanol1,method = "euclidean",diag = T)
# Ethanol1_dis <- Ethanol1_dis %>% as.matrix %>% data.frame(check.names = F,stringsAsFactors = F)
# diag(Ethanol1_dis) <- NA ;Ethanol1_dis[upper.tri(Ethanol1_dis)] <- NA
# Ethanol1_dis <- Ethanol1_dis %>% rownames_to_column("sample1")
# Ethanol1_dis <- Ethanol1_dis %>% pivot_longer(cols = -1,names_to="sample2",values_to="Ethanol") %>% filter(!is.na(Ethanol))
# Ethanol1_dis <- left_join(Ethanol1_dis,group,by=c("sample1"="sampleid")) %>% left_join(group,by=c("sample2"="sampleid"),suffix=c("1","2"))
# Ethanol1_dis <- filter(Ethanol1_dis,StarterType1==StarterType2)
# Ethanol1_dist <- Ethanol1_dis[1:3] %>%
#   pivot_wider(names_from = "sample2",values_from="Ethanol") %>%
#   column_to_rownames("sample1") %>% as.dist
# 
# bnti_Ethanol1 <- left_join(bac_bnti1_dis[1:3],Ethanol1_dis)
# summary(lm(betanti~Ethanol,data = bnti_Ethanol1))
# ggplot(bnti_Ethanol1,aes(Ethanol,betanti))+geom_point()+geom_smooth(method = "lm")
# 
# mantel(Ethanol1_dist,bac_bnti1_dist,na.rm = T)
# 
# #  phase 3
# 
# bac_bnti3_dis <- bac_bnti_dist[c(22:33,55:66),c(22:33,55:66)]
# bac_bnti3_dis <- bac_bnti3_dis %>% rownames_to_column("sample1")
# bac_bnti3_dis <- bac_bnti3_dis %>% pivot_longer(cols = -1,names_to="sample2",values_to="betanti") %>% filter(!is.na(betanti))
# bac_bnti3_dis <- left_join(bac_bnti3_dis,group,by=c("sample1"="sampleid")) %>% left_join(group,by=c("sample2"="sampleid"),suffix=c("1","2"))
# bac_bnti3_dis <- filter(bac_bnti3_dis,StarterType1==StarterType2)
# bac_bnti3_dist <- bac_bnti3_dis[1:3] %>%
#   pivot_wider(names_from = "sample2",values_from="betanti") %>%
#   column_to_rownames("sample1") %>% as.dist
# 
# # Moisture
# Moisture3 <- env3 %>% select(sampleid,Moisture) %>% column_to_rownames("sampleid")
# Moisture3_dis <- vegdist(Moisture3,method = "euclidean",diag = T)
# 
# Moisture3_dis <- Moisture3_dis %>% as.matrix %>% data.frame(check.names = F,stringsAsFactors = F)
# diag(Moisture3_dis) <- NA ;Moisture3_dis[upper.tri(Moisture3_dis)] <- NA
# Moisture3_dis <- Moisture3_dis %>% rownames_to_column("sample1")
# Moisture3_dis <- Moisture3_dis %>% pivot_longer(cols = -1,names_to="sample2",values_to="Moisture") %>% filter(!is.na(Moisture))
# Moisture3_dis <- left_join(Moisture3_dis,group,by=c("sample1"="sampleid")) %>% left_join(group,by=c("sample2"="sampleid"),suffix=c("1","2"))
# Moisture3_dis <- filter(Moisture3_dis,StarterType1==StarterType2)
# Moisture3_dist <- Moisture3_dis[1:3] %>%
#   pivot_wider(names_from = "sample2",values_from="Moisture") %>%
#   column_to_rownames("sample1") %>% as.dist
# 
# bnti_Moisture3 <- left_join(bac_bnti3_dis[1:3],Moisture3_dis)
# summary(lm(betanti~Moisture,data = bnti_Moisture3))
# ggplot(bnti_Moisture3,aes(Moisture,betanti))+geom_point()+geom_smooth(method = "lm")
# 
# mantel(Moisture3_dist,bac_bnti3_dist,na.rm=T)
# 
# # Reducing Sugar
# ReducingSugar3 <- env3 %>% select(sampleid,ReducingSugar) %>% column_to_rownames("sampleid")
# ReducingSugar3_dis <- vegdist(ReducingSugar3,method = "euclidean",diag = T)
# 
# ReducingSugar3_dis <- ReducingSugar3_dis %>% as.matrix %>% data.frame(check.names = F,stringsAsFactors = F)
# diag(ReducingSugar3_dis) <- NA ;ReducingSugar3_dis[upper.tri(ReducingSugar3_dis)] <- NA
# ReducingSugar3_dis <- ReducingSugar3_dis %>% rownames_to_column("sample1")
# ReducingSugar3_dis <- ReducingSugar3_dis %>% pivot_longer(cols = -1,names_to="sample2",values_to="ReducingSugar") %>% filter(!is.na(ReducingSugar))
# ReducingSugar3_dis <- left_join(ReducingSugar3_dis,group,by=c("sample1"="sampleid")) %>% left_join(group,by=c("sample2"="sampleid"),suffix=c("1","2"))
# ReducingSugar3_dis <- filter(ReducingSugar3_dis,StarterType1==StarterType2)
# ReducingSugar3_dist <- ReducingSugar3_dis[1:3] %>%
#   pivot_wider(names_from = "sample2",values_from="ReducingSugar") %>%
#   column_to_rownames("sample1") %>% as.dist
# 
# 
# bnti_ReducingSugar3 <- left_join(bac_bnti3_dis[1:3],ReducingSugar3_dis)
# summary(lm(betanti~ReducingSugar,data = bnti_ReducingSugar3))
# ggplot(bnti_ReducingSugar3,aes(ReducingSugar,betanti))+geom_point()+geom_smooth(method = "lm")
# 
# mantel(ReducingSugar3_dist,bac_bnti3_dist,na.rm=T)
# 
# # Acidity
# 
# Acidity3 <- env3 %>% select(sampleid,Acidity) %>% column_to_rownames("sampleid")
# Acidity3_dis <- vegdist(Acidity3,method = "euclidean",diag = T)
# 
# Acidity3_dis <- Acidity3_dis %>% as.matrix %>% data.frame(check.names = F,stringsAsFactors = F)
# diag(Acidity3_dis) <- NA ;Acidity3_dis[upper.tri(Acidity3_dis)] <- NA
# Acidity3_dis <- Acidity3_dis %>% rownames_to_column("sample1")
# Acidity3_dis <- Acidity3_dis %>% pivot_longer(cols = -1,names_to="sample2",values_to="Acidity") %>% filter(!is.na(Acidity))
# Acidity3_dis <- left_join(Acidity3_dis,group,by=c("sample1"="sampleid")) %>% left_join(group,by=c("sample2"="sampleid"),suffix=c("1","2"))
# Acidity3_dis <- filter(Acidity3_dis,StarterType1==StarterType2)
# Acidity3_dist <- Acidity3_dis[1:3] %>%
#   pivot_wider(names_from = "sample2",values_from="Acidity") %>%
#   column_to_rownames("sample1") %>% as.dist
# 
# bnti_Acidity3 <- left_join(bac_bnti3_dis[1:3],Acidity3_dis)
# summary(lm(betanti~Acidity,data = bnti_Acidity3))
# ggplot(bnti_Acidity3,aes(Acidity,betanti))+geom_point()+geom_smooth(method = "lm")
# 
# mantel(Acidity3_dist,bac_bnti3_dist,na.rm=T)
# 
# # Temperature
# 
# Temperature3 <- env3 %>% select(sampleid,Temperature) %>% column_to_rownames("sampleid")
# Temperature3_dis <- vegdist(Temperature3,method = "euclidean",diag = T)
# Temperature3_dis <- Temperature3_dis %>% as.matrix %>% data.frame(check.names = F,stringsAsFactors = F)
# diag(Temperature3_dis) <- NA ;Temperature3_dis[upper.tri(Temperature3_dis)] <- NA
# Temperature3_dis <- Temperature3_dis %>% rownames_to_column("sample1")
# Temperature3_dis <- Temperature3_dis %>% pivot_longer(cols = -1,names_to="sample2",values_to="Temperature") %>% filter(!is.na(Temperature))
# Temperature3_dis <- left_join(Temperature3_dis,group,by=c("sample1"="sampleid")) %>% left_join(group,by=c("sample2"="sampleid"),suffix=c("1","2"))
# Temperature3_dis <- filter(Temperature3_dis,StarterType1==StarterType2)
# Temperature3_dist <- Temperature3_dis[1:3] %>%
#   pivot_wider(names_from = "sample2",values_from="Temperature") %>%
#   column_to_rownames("sample1") %>% as.dist
# 
# bnti_Temperature3 <- left_join(bac_bnti3_dis[1:3],Temperature3_dis)
# summary(lm(betanti~Temperature,data = bnti_Temperature3))
# ggplot(bnti_Temperature3,aes(Temperature,betanti))+geom_point()+geom_smooth(method = "lm")
# 
# mantel(Temperature3_dist,bac_bnti3_dist,na.rm=T)
# 
# # pH
# 
# pH3 <- env3 %>% select(sampleid,pH) %>% column_to_rownames("sampleid")
# pH3_dis <- vegdist(pH3,method = "euclidean",diag = T)
# pH3_dis <- pH3_dis %>% as.matrix %>% data.frame(check.names = F,stringsAsFactors = F)
# diag(pH3_dis) <- NA ;pH3_dis[upper.tri(pH3_dis)] <- NA
# pH3_dis <- pH3_dis %>% rownames_to_column("sample1")
# pH3_dis <- pH3_dis %>% pivot_longer(cols = -1,names_to="sample2",values_to="pH") %>% filter(!is.na(pH))
# pH3_dis <- left_join(pH3_dis,group,by=c("sample1"="sampleid")) %>% left_join(group,by=c("sample2"="sampleid"),suffix=c("1","2"))
# pH3_dis <- filter(pH3_dis,StarterType1==StarterType2)
# pH3_dist <- pH3_dis[1:3] %>%
#   pivot_wider(names_from = "sample2",values_from="pH") %>%
#   column_to_rownames("sample1") %>% as.dist
# 
# bnti_pH3 <- left_join(bac_bnti3_dis[1:3],pH3_dis)
# summary(lm(betanti~pH,data = bnti_pH3))
# ggplot(bnti_pH3,aes(pH,betanti))+geom_point()+geom_smooth(method = "lm")
# 
# mantel(pH3_dist,bac_bnti3_dist,na.rm = T)
# 
# # Ethanol
# Ethanol3 <- env3 %>% select(sampleid,Ethanol) %>% column_to_rownames("sampleid")
# Ethanol3_dis <- vegdist(Ethanol3,method = "euclidean",diag = T)
# Ethanol3_dis <- Ethanol3_dis %>% as.matrix %>% data.frame(check.names = F,stringsAsFactors = F)
# diag(Ethanol3_dis) <- NA ;Ethanol3_dis[upper.tri(Ethanol3_dis)] <- NA
# Ethanol3_dis <- Ethanol3_dis %>% rownames_to_column("sample1")
# Ethanol3_dis <- Ethanol3_dis %>% pivot_longer(cols = -1,names_to="sample2",values_to="Ethanol") %>% filter(!is.na(Ethanol))
# Ethanol3_dis <- left_join(Ethanol3_dis,group,by=c("sample1"="sampleid")) %>% left_join(group,by=c("sample2"="sampleid"),suffix=c("1","2"))
# Ethanol3_dis <- filter(Ethanol3_dis,StarterType1==StarterType2)
# Ethanol3_dist <- Ethanol3_dis[1:3] %>%
#   pivot_wider(names_from = "sample2",values_from="Ethanol") %>%
#   column_to_rownames("sample1") %>% as.dist
# 
# bnti_Ethanol3 <- left_join(bac_bnti3_dis[1:3],Ethanol3_dis)
# summary(lm(betanti~Ethanol,data = bnti_Ethanol3))
# ggplot(bnti_Ethanol3,aes(Ethanol,betanti))+geom_point()+geom_smooth(method = "lm")
# 
# mantel(Ethanol3_dist,bac_bnti3_dist,na.rm = T)


# RCBC --------------------------------------------------------------------



# 16s

bac_betanti <- read.delim("16s/betanti/bacteria_betaNTI.tsv",row.names = 1) %>% as.matrix()

bac_rcbc <- read.delim("16s/betanti/bacteria_RCBC.tsv",row.names = 1)
diag(bac_rcbc) <- NA;bac_rcbc[upper.tri(bac_rcbc)] <- NA
bac_rcbc <- bac_rcbc %>% as.matrix()

bac_rcbc[abs(bac_betanti)>=2] <- NA

bac_rcbc <- bac_rcbc %>% data.frame %>% 
  rownames_to_column("sample1") %>% 
  pivot_longer(cols=-1,names_to="sample2",values_to="RCbray") %>% 
  filter(!is.na(RCbray))

bac_rcbc <- left_join(bac_rcbc,group,by=c("sample1"="sampleid")) %>% 
  left_join(group,by=c("sample2"="sampleid"),suffix=c("1","2"))

bac_rcbc <- filter(bac_rcbc,SampleTime1==SampleTime2)

# whole

  bac_rcbc %>% ggplot(aes(SampleTime1,RCbray))+
  geom_hline(yintercept = c(-0.95,0.95),linetype="dashed")+
  geom_boxplot(outlier.colour = NA)+geom_jitter(aes(color=StarterType1))+
  scale_y_continuous(breaks = c(-1,-0.95,-0.5,0,0.5,0.95,1.0))+
  scale_x_discrete(limits=c("SAC0h","SAC6h","SAC12h","SAC18h","SAC24h",
                            "FER0d","FER3d","FER6d","FER9d","FER12d","FER15d"))+
  theme_bw(base_size = 8)+theme(panel.grid = element_blank(),
                   axis.title.x = element_blank(),
                   axis.text.x = element_text(angle = 90,vjust = 0.5))+
    scale_color_manual(values=c("PureStarter"="#00BFC4","TraditionalStarter"="#F8766D"))+
    labs(color="Batch")
  
  
  
  its_rcbc %>% ggplot(aes(SampleTime1,RCbray))+
    geom_hline(yintercept = c(-0.95,0.95),linetype="dashed")+
    geom_boxplot()+geom_jitter(aes(color=StarterType1))+
    scale_y_continuous(breaks = c(-1,-0.95,-0.5,0,0.5,0.95,1.0))+
    scale_x_discrete(limits=c("SAC0h","SAC6h","SAC12h","SAC18h","SAC24h",
                              "FER0d","FER3d","FER6d","FER9d","FER12d","FER15d"))+
    theme_bw()+theme(panel.grid = element_blank(),
                     axis.title.x = element_blank(),
                     axis.text.x = element_text(angle = 90,vjust = 0.5))+
    scale_color_manual(values=c("PureStarter"="#00BFC4","TraditionalStarter"="#F8766D"))+
    labs(color="Batch")
  
  
  
  
  
# pure
bac_rcbc %>% filter(StarterType1==StarterType2,StarterType1=="PureStarter") %>% 
  ggplot(aes(SampleTime1,RCbray))+
  geom_hline(yintercept = c(-0.95,0.95),linetype="dashed")+
  geom_boxplot(outlier.colour = NA)+geom_jitter()+
  scale_y_continuous(breaks = c(-1,-0.95,-0.5,0,0.5,0.95,1.0))+
  scale_x_discrete(limits=c("SAC0h","SAC6h","SAC12h","SAC18h","SAC24h",
                            "FER0d","FER3d","FER6d","FER9d","FER12d","FER15d"))+
  theme_bw()+theme(panel.grid = element_blank(),
                   axis.title.x = element_blank())
# traditional

bac_rcbc %>% filter(StarterType1==StarterType2,StarterType1=="TraditionalStarter") %>% 
  ggplot(aes(SampleTime1,RCbray))+
  geom_hline(yintercept = c(-0.95,0.95),linetype="dashed")+
  geom_boxplot(outlier.colour = NA)+geom_jitter()+
  scale_y_continuous(breaks = c(-1,-0.95,-0.5,0,0.5,0.95,1.0))+
  scale_x_discrete(limits=c("SAC0h","SAC6h","SAC12h","SAC18h","SAC24h",
                            "FER0d","FER3d","FER6d","FER9d","FER12d","FER15d"))+
  theme_bw()+theme(panel.grid = element_blank(),
                   axis.title.x = element_blank())


# its

its_betanti <- read.delim("its/betanti/fungi_betaNTI.tsv",row.names = 1) %>% as.matrix()

its_rcbc <- read.delim("its/betanti/fungi_RCBC.tsv",row.names = 1)
diag(its_rcbc) <- NA;its_rcbc[upper.tri(its_rcbc)] <- NA
its_rcbc <- its_rcbc %>% as.matrix()

its_rcbc[abs(its_betanti)>=2] <- NA

its_rcbc <- its_rcbc %>% data.frame %>% 
  rownames_to_column("sample1") %>% 
  pivot_longer(cols=-1,names_to="sample2",values_to="RCbray") %>% 
  filter(!is.na(RCbray))

its_rcbc <- left_join(its_rcbc,group,by=c("sample1"="sampleid")) %>% 
  left_join(group,by=c("sample2"="sampleid"),suffix=c("1","2"))

its_rcbc <- filter(its_rcbc,SampleTime1==SampleTime2)

# whole

its_rcbc %>% ggplot(aes(SampleTime1,RCbray))+
  geom_hline(yintercept = c(-0.95,0.95),linetype="dashed")+
  geom_boxplot()+geom_jitter(aes(color=StarterType1))+
  scale_y_continuous(breaks = c(-1,-0.95,-0.5,0,0.5,0.95,1.0))+
  scale_x_discrete(limits=c("SAC0h","SAC6h","SAC12h","SAC18h","SAC24h",
                            "FER0d","FER3d","FER6d","FER9d","FER12d","FER15d"))+
  theme_bw()+theme(panel.grid = element_blank(),
                   axis.title.x = element_blank(),
                   axis.text.x = element_text(angle = 90,vjust = 0.5))+
  scale_color_manual(values=c("PureStarter"="#00BFC4","TraditionalStarter"="#F8766D"))+
  labs(color="Batch")

# pure
its_rcbc %>% filter(StarterType1==StarterType2,StarterType1=="PureStarter") %>% 
  ggplot(aes(SampleTime1,RCbray))+
  geom_hline(yintercept = c(-0.95,0.95),linetype="dashed")+
  geom_boxplot()+geom_point()+
  scale_y_continuous(breaks = c(-1,-0.95,-0.5,0,0.5,0.95,1.0))+
  scale_x_discrete(limits=c("SAC0h","SAC6h","SAC12h","SAC18h","SAC24h",
                            "FER0d","FER3d","FER6d","FER9d","FER12d","FER15d"))+
  theme_bw()+theme(panel.grid = element_blank(),
                   axis.title.x = element_blank())
# traditional

its_rcbc %>% filter(StarterType1==StarterType2,StarterType1=="TraditionalStarter") %>% 
  ggplot(aes(SampleTime1,RCbray))+
  geom_hline(yintercept = c(-0.95,0.95),linetype="dashed")+
  geom_boxplot()+geom_point()+
  scale_y_continuous(breaks = c(-1,-0.95,-0.5,0,0.5,0.95,1.0))+
  scale_x_discrete(limits=c("SAC0h","SAC6h","SAC12h","SAC18h","SAC24h",
                            "FER0d","FER3d","FER6d","FER9d","FER12d","FER15d"))+
  theme_bw()+theme(panel.grid = element_blank(),
                   axis.title.x = element_blank())

