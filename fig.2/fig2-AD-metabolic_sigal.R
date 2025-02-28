library(tidyverse)



metabolic <- read.csv("sigal_data.csv")

metabolic1 <- metabolic %>% 
  group_by(Inoculate,time) %>% 
  mutate(Mois_mean=mean(Moisture),
         Mois_sd=sd(Moisture),
         pH_mean=mean(pH),
         pH_sd=sd(pH),
         acidity_mean=mean(Acidity),
         acidity_sd=sd(Acidity),
         Redu_mean=mean(Reducing_sugar),
         Redu_sd=sd(Reducing_sugar),
         .before=2)
metabolic1$Inoculate <- factor(metabolic1$Inoculate,levels = unique(metabolic1$Inoculate), ordered = T)
metabolic1$time <- factor(metabolic1$time,levels = unique(metabolic1$time))

colnames(metabolic1)
# [1] "Sample"         "Mois_mean"      "Mois_sd"        "pH_mean"        "pH_sd"          "acidity_mean"   "acidity_sd"    
# [8] "Redu_mean"      "Redu_sd"        "Inoculate"      "time"           "Moisture"       "pH"             "Acidity"       
# [15] "Reducing_sugar"


#满足假设，单因素方差分析，详情使用?aov查看帮助
#Moisture
fit <- aov(Moisture~time, data = metabolic1[1:15,])
summary(fit)#no signal
fit <- aov(Moisture~time, data = metabolic1[16:30,])
summary(fit)#no signal
fit <- aov(Moisture~time, data = metabolic1[31:45,])
summary(fit)#no signal
fit <- aov(Moisture~time, data = metabolic1[46:60,])
summary(fit)#no signal
fit <- aov(Moisture~time, data = metabolic1[61:75,])
summary(fit)#no signal
fit <- aov(Moisture~time, data = metabolic1[76:90,])
summary(fit)#signal
#Tukey HSD 检验http://127.0.0.1:38985/graphics/plot_zoom_png?width=1200&height=900
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk
library(multcomp)#multcomp包中提供了更直观的方法，展示Tukey检验的结果。
tuk <- glht(fit, alternative = 'two.sided', linfct = mcp(site = 'Tukey'))
plot(cld(tuk, level = 0.05, decreasing = TRUE))


#pH
fit <- aov(pH~time, data = metabolic1[1:15,])
summary(fit)#signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk
fit <- aov(pH~time, data = metabolic1[16:30,])
summary(fit)#signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk
fit <- aov(pH~time, data = metabolic1[31:45,])
summary(fit)#signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk
fit <- aov(pH~time, data = metabolic1[46:60,])
summary(fit)#no signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk
fit <- aov(pH~time, data = metabolic1[61:75,])
summary(fit)#signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk
fit <- aov(pH~time, data = metabolic1[76:90,])
summary(fit)#signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk

#Acidity
fit <- aov(Acidity~time, data = metabolic1[1:15,])
summary(fit)#signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk
fit <- aov(Acidity~time, data = metabolic1[16:30,])
summary(fit)#signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk
fit <- aov(Acidity~time, data = metabolic1[31:45,])
summary(fit)#signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk
fit <- aov(Acidity~time, data = metabolic1[46:60,])
summary(fit)#no signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk
fit <- aov(Acidity~time, data = metabolic1[61:75,])
summary(fit)#signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk
fit <- aov(Acidity~time, data = metabolic1[76:90,])
summary(fit)#signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk


#Reducing_sugar
fit <- aov(Reducing_sugar~time, data = metabolic1[1:15,])
summary(fit)#signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk
fit <- aov(Reducing_sugar~time, data = metabolic1[16:30,])
summary(fit)#signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk
fit <- aov(Reducing_sugar~time, data = metabolic1[31:45,])
summary(fit)#signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk
fit <- aov(Reducing_sugar~time, data = metabolic1[46:60,])
summary(fit)#no signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk
fit <- aov(Reducing_sugar~time, data = metabolic1[61:75,])
summary(fit)#signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk
fit <- aov(Reducing_sugar~time, data = metabolic1[76:90,])
summary(fit)#signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk


#满足假设，单因素方差分析，详情使用?aov查看帮助
fit <- aov(Acidity~time, data = metabolic1)
summary(fit)
#Tukey HSD 检验
tuk <- TukeyHSD(fit, conf.level = 0.95) 
plot(tuk)#没有碰到虚线的为显著关系
tuk
#满足假设，单因素方差分析，详情使用?aov查看帮助
fit <- aov(Reducing_sugar~time, data = metabolic1)


#添加水分多重比较的结果
metabolic1 <- metabolic1[-(16:18)]
metabolic1$sig_mois <- c(rep("a",78),rep("ab",3),rep("b",3),rep("ab",3),rep("a",3))
metabolic1$sig_ph <- c(rep("a",3),rep("c",3),rep("bc",3),rep("bc",3),rep("b",3),
            rep("a",3),rep("d",3),rep("bc",3),rep("cd",3),rep("b",3),
            rep("a",3),rep("d",3),rep("c",3),rep("cd",3),rep("b",3),
            rep("a",3),rep("d",3),rep("d",3),rep("c",3),rep("b",3),
            rep("a",3),rep("cd",3),rep("d",3),rep("c",3),rep("b",3),
            rep("a",3),rep("c",3),rep("c",3),rep("a",3),rep("b",3)
            )
metabolic1$sig_aci <- c(rep("d",3),rep("c",3),rep("bc",3),rep("a",3),rep("b",3),
                       rep("f",3),rep("d",3),rep("b",3),rep("a",3),rep("c",3),
                       rep("f",3),rep("d",3),rep("b",3),rep("a",3),rep("c",3),
                       rep("d",3),rep("c",3),rep("ab",3),rep("a",3),rep("bc",3),
                       rep("cd",3),rep("c",3),rep("ab",3),rep("a",3),rep("bc",3),
                       rep("cd",3),rep("c",3),rep("ab",3),rep("a",3),rep("bc",3)
)
metabolic1$sig_sugare <- c(rep("a",3),rep("c",3),rep("b",3),rep("bc",3),rep("bc",3),
                          rep("a",3),rep("b",3),rep("ab",3),rep("b",3),rep("b",3),
                          rep("a",3),rep("b",3),rep("a",3),rep("b",3),rep("b",3),
                          rep("a",3),rep("c",3),rep("a",3),rep("bc",3),rep("ab",3),
                          rep("a",3),rep("c",3),rep("b",3),rep("c",3),rep("c",3),
                          rep("a",3),rep("c",3),rep("b",3),rep("bc",3),rep("bc",3)
                          )
write.csv(metabolic1,"result/metabolic1.csv")

library(ggsci)
cols <- pal_npg()(8)
metabolic1 <- metabolic1 %>% 
  mutate(Inoculate = recode(Inoculate,"LAB" = "Lactobacillus", "WESS" = "Weissella",
                            "Pichi" = "Pichia", "Sac" = "Saccharomyces", "Canda" = "Candida"))

p1 <- ggplot(metabolic1,aes(time,Mois_mean,color=Inoculate))+
  geom_point()+
  geom_line(aes(group=Inoculate))+
  geom_errorbar(aes(ymin=Mois_mean-Mois_sd,ymax=Mois_mean+Mois_sd),width = .2)+
  scale_color_npg()+
  scale_y_continuous(limits = c(30,50),n.breaks = 6)+
  theme_bw(base_size = 8)+theme(panel.grid = element_blank(),
                                legend.position = "right",
                                legend.text = element_text(color="black",size=8,face="italic"), # 定义图例文本
                                legend.justification =c(1,0.5) )+
  facet_wrap(Inoculate~.)+
  geom_text(aes(label = sig_mois, y = 50),color = "black") +
  xlab(NULL)+ylab("Moisture content %")+
  theme(axis.title.y = element_text(size = 14),
        strip.text = element_text(size = rel(1.5),face = "italic"),
        strip.background = element_rect(fill = "lightblue", colour ="black",
        size = 1)
        )
p1
ggsave(p1,device = "png",file="result/Moisture.png") 
ggsave(p1,device = "pdf",file="result/Moisture content.pdf") 

p2 <- ggplot(metabolic1,aes(time,pH_mean,color=Inoculate))+
  geom_point()+
  geom_line(aes(group=Inoculate))+
  geom_errorbar(aes(ymin=pH_mean-pH_sd,ymax=pH_mean+pH_sd),width = .2)+
  scale_color_npg()+
  scale_y_continuous(limits = c(3.5,4.65),n.breaks = 6)+
  theme_bw(base_size = 8)+theme(panel.grid = element_blank(),
                                legend.position = "right",
                                legend.text = element_text(color="black",size=8,face="italic"), # 定义图例文本
                                legend.justification =c(1,0.5) )+
  facet_wrap(Inoculate~.)+
  geom_text(aes(label = sig_ph, y = 4.6),color = "black") +
  xlab(NULL)+ylab("pH")+
  theme(axis.title.y = element_text(size = 14),
        strip.text = element_text(size = rel(1.5),face = "italic"),
        strip.background = element_rect(fill = "lightblue", colour ="black",
                                        size = 1)
  )
p2
ggsave(p2,device = "png",file="result/pH.png") 
ggsave(p2,device = "pdf",file="result/pH.pdf") 

p3 <- ggplot(metabolic1,aes(time,acidity_mean,color=Inoculate))+
  geom_point()+
  geom_line(aes(group=Inoculate))+
  geom_errorbar(aes(ymin=acidity_mean-acidity_sd,ymax=acidity_mean+acidity_sd),width = .2)+
  scale_color_npg()+
  scale_y_continuous(limits = c(0,10),n.breaks = 6)+
  theme_bw(base_size = 8)+theme(panel.grid = element_blank(),
                                legend.position = "right",
                                legend.text = element_text(color="black",size=8,face="italic"), # 定义图例文本
                                legend.justification =c(1,0.5) )+
  facet_wrap(Inoculate~.)+
  geom_text(aes(label = sig_aci, y = 9),color = "black") +
  xlab(NULL)+ylab("Acidity (°T)")+
  theme(axis.title.y = element_text(size = 14),
        strip.text = element_text(size = rel(1.5),face = "italic"),
        strip.background = element_rect(fill = "lightblue", colour ="black",
                                        size = 1)
  )
p3
ggsave(p3,device = "png",file="result/acidity.png") 
ggsave(p3,device = "pdf",file="result/acidity.pdf")

p4 <- ggplot(metabolic1,aes(time,Redu_mean,color=Inoculate))+
  geom_point()+
  geom_line(aes(group=Inoculate))+
  geom_errorbar(aes(ymin=Redu_mean-Redu_sd,ymax=Redu_mean+Redu_sd),width = .2)+
  scale_color_npg()+
  scale_y_continuous(limits = c(0,10),n.breaks = 6)+
  theme_bw(base_size = 8)+theme(panel.grid = element_blank(),
                                legend.position = "right",
                                legend.text = element_text(color="black",size=8,face="italic"), # 定义图例文本
                                legend.justification =c(1,0.5) )+
  facet_wrap(Inoculate~.)+
  geom_text(aes(label = sig_sugare, y = 9),color = "black") +
  xlab(NULL)+ylab("Reducing sugar(g/kg)")+
  theme(axis.title.y = element_text(size = 14),
        strip.text = element_text(size = rel(1.5),face = "italic"),
        strip.background = element_rect(fill = "lightblue", colour ="black",
                                        size = 1)
  )
p4
ggsave(p4,device = "png",file="result/Reducing sugar.png") 
ggsave(p4,device = "pdf",file="result/Reducing sugar.pdf")

library(patchwork)
p1 <- p1 + 
  theme(axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = rel(1.5),face = "italic"),
        strip.background = element_rect(fill = "lightblue", colour ="black",
                                        size = 1)
  )
p2 <- p2 + 
  theme(axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = rel(1.5),face = "italic"),
        strip.background = element_rect(fill = "lightblue", colour ="black",
                                        size = 1)
  )
p3 <- p3 + 
  theme(axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = rel(1.5),face = "italic"),
        strip.background = element_rect(fill = "lightblue", colour ="black",
                                        size = 1)
  )
p4 <- p4 + 
  theme(axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = rel(1.5),face = "italic"),
        strip.background = element_rect(fill = "lightblue", colour ="black",
                                        size = 1)
  )
b1 <- (p1+p2+p3+p4)+plot_layout(guides = 'collect')&
  theme(legend.position='bottom')
b1
ggsave(b1 ,file="result/lihua.pdf",width = 24 , height = 14,units = "cm")
ggsave(b1 ,file="result/lihua.png",width = 24 , height = 14,units = "cm",dpi = 300)


#----------------------------------------test
library(dplyr)
library(ggplot2)
library(ggsci)
mpg
metabolic1 %>% pivot_longer(cols = !Sample,names_to = "Samples",
                    values_to = "number") 

acid <- read.csv("yexiang_data.csv")
colnames(acid)
# "Sample"      "Inoculate"   "time"        "Lactic_acid" "Acetic_acid" "Ethanol"     "Glucose" 

acid1 <- acid %>% 
  group_by(Inoculate,time) %>% 
  mutate(lac_mean=mean(Lactic_acid),
         lac_sd=sd(Lactic_acid),
         ace_mean=mean(Acetic_acid),
         ace_sd=sd(Acetic_acid),
         eth_mean=mean(Ethanol),
         eth_sd=sd(Ethanol),
         glu_mean=mean(Glucose),
         glu_sd=sd(Glucose),
         .before=2)
acid1$Inoculate <- factor(acid1$Inoculate,levels = unique(acid1$Inoculate), ordered = T)
acid1$time <- factor(acid1$time,levels = unique(acid1$time))

colnames(acid1)
# [1] "Sample"      "lac_mean"    "lac_sd"      "ace_mean"    "ace_sd"      "eth_mean"    "eth_sd"      "glu_mean"    "glu_sd"     
# [10] "Inoculate"   "time"        "Lactic_acid" "Acetic_acid" "Ethanol"     "Glucose"    

#Lactic_acid
fit <- aov(Lactic_acid~time, data = acid1[1:15,])
summary(fit)#signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk
fit <- aov(Lactic_acid~time, data = acid1[16:30,])
summary(fit)#signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk
fit <- aov(Lactic_acid~time, data = acid1[31:45,])
summary(fit)#signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk
fit <- aov(Lactic_acid~time, data = acid1[46:60,])
summary(fit)#signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk
fit <- aov(Lactic_acid~time, data = acid1[61:75,])
summary(fit)#signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk
fit <- aov(Lactic_acid~time, data = acid1[76:90,])
summary(fit)#signal
#Tukey HSD 检验
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk
library(multcomp)#multcomp包中提供了更直观的方法，展示Tukey检验的结果。
tuk <- glht(fit, alternative = 'two.sided', linfct = mcp(site = 'Tukey'))
plot(cld(tuk, level = 0.05, decreasing = TRUE))


#Acetic_acid
fit <- aov(Acetic_acid~time, data = acid1[1:15,])
summary(fit)#signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk
fit <- aov(Acetic_acid~time, data = acid1[16:30,])
summary(fit)#signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk
fit <- aov(Acetic_acid~time, data = acid1[31:45,])
summary(fit)#signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk
fit <- aov(Acetic_acid~time, data = acid1[46:60,])
summary(fit)#signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk
fit <- aov(Acetic_acid~time, data = acid1[61:75,])
summary(fit)#signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk
fit <- aov(Acetic_acid~time, data = acid1[76:90,])
summary(fit)#signal
#Tukey HSD 检验
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk


#Ethanol
fit <- aov(Ethanol~time, data = acid1[1:15,])
summary(fit)#no signal

fit <- aov(Ethanol~time, data = acid1[16:30,])
summary(fit)#no signal

fit <- aov(Ethanol~time, data = acid1[31:45,])
summary(fit)#signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk
fit <- aov(Ethanol~time, data = acid1[46:60,])
summary(fit)#no signal

fit <- aov(Ethanol~time, data = acid1[61:75,])
summary(fit)#signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk
fit <- aov(Ethanol~time, data = acid1[76:90,])
summary(fit)#signal
#Tukey HSD 检验
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk


#Glucose
fit <- aov(Glucose~time, data = acid1[1:15,])
summary(fit)#signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk
fit <- aov(Glucose~time, data = acid1[16:30,])
summary(fit)#signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk
fit <- aov(Glucose~time, data = acid1[31:45,])
summary(fit)#signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk
fit <- aov(Glucose~time, data = acid1[46:60,])
summary(fit)#no signal

fit <- aov(Glucose~time, data = acid1[61:75,])
summary(fit)#signal
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk
fit <- aov(Glucose~time, data = acid1[76:90,])
summary(fit)#signal
#Tukey HSD 检验
tuk <- TukeyHSD(fit, conf.level = 0.95)
plot(tuk)#没有碰到虚线的为显著关系
tuk


#添加水分多重比较的结果
acid1$sig_lac <- c(rep("d",3),rep("cd",3),rep("c",3),rep("a",3),rep("b",3),
                   rep("c",3),rep("bc",3),rep("a",3),rep("a",3),rep("ab",3),
                   rep("c",3),rep("bc",3),rep("a",3),rep("b",3),rep("b",3),
                   rep("b",3),rep("ab",3),rep("ab",3),rep("a",3),rep("a",3),
                   rep("c",3),rep("b",3),rep("ab",3),rep("a",3),rep("ab",3),
                   rep("d",3),rep("c",3),rep("a",3),rep("b",3),rep("bc",3)
)

acid1$sig_ace <- c(rep("b",3),rep("b",3),rep("b",3),rep("a",3),rep("a",3),
                   rep("b",3),rep("b",3),rep("b",3),rep("ab",3),rep("a",3),
                   rep("b",3),rep("b",3),rep("b",3),rep("ab",3),rep("a",3),
                   rep("b",3),rep("b",3),rep("b",3),rep("ab",3),rep("a",3),
                   rep("ab",3),rep("b",3),rep("ab",3),rep("ab",3),rep("a",3),
                   rep("c",3),rep("c",3),rep("bc",3),rep("ab",3),rep("a",3)
)

acid1$sig_eth <- c(rep("a",15),
                   rep("a",15),
                   rep("b",3),rep("ab",3),rep("a",3),rep("b",3),rep("b",3),
                   rep("a",15),
                   rep("a",3),rep("b",3),rep("b",3),rep("b",3),rep("b",3),
                   rep("d",3),rep("b",3),rep("a",3),rep("b",3),rep("c",3)
)

acid1$sig_glu <- c(rep("a",3),rep("b",3),rep("b",3),rep("b",3),rep("b",3),
                   rep("ab",3),rep("a",3),rep("b",3),rep("b",3),rep("ab",3),
                   rep("ab",3),rep("b",3),rep("a",3),rep("ab",3),rep("ab",3),
                   rep("a",15),
                   rep("a",3),rep("b",3),rep("b",3),rep("b",3),rep("b",3),
                   rep("a",3),rep("b",3),rep("bc",3),rep("bc",3),rep("c",3)
)
write.csv(acid1,"result/acid1.csv")

acid1 <- read.csv("result/acid1.csv", row.names = 1)

acid1 <- acid1 %>% 
  mutate(Inoculate = recode(Inoculate,"LAB" = "Lactobacillus", "WESS" = "Weissella",
                            "Pichi" = "Pichia", "Sac" = "Saccharomyces", "Canda" = "Candida"))
acid1$time <- factor(acid1$time,levels = unique(acid1$time))

q1 <- ggplot(acid1,aes(time,lac_mean,color=Inoculate))+
  geom_point()+
  geom_line(aes(group=Inoculate))+
  geom_errorbar(aes(ymin=lac_mean-lac_sd,ymax=lac_mean+lac_sd),width = .2)+
  scale_color_npg()+
  scale_y_continuous(limits = c(0,10),n.breaks = 6)+
  theme_bw(base_size = 8)+theme(panel.grid = element_blank(),
                                legend.position = "right",
                                legend.text = element_text(color="black",size=8,face="italic"), # 定义图例文本
                                legend.justification =c(1,0.5) )+
  facet_wrap(Inoculate~.)+
  geom_text(aes(label = sig_lac, y = 9),color = "black") +
  xlab(NULL)+ylab("Lactic acid content (mg/g)")+
  theme(axis.title.y = element_text(size = 14),
        strip.text = element_text(size = rel(1.5),face = "italic"),
        strip.background = element_rect(fill = "lightblue", colour ="black",
                                        size = 1)
  )
q1
ggsave(q1,device = "png",file="result/lac.png") 
ggsave(q1,device = "pdf",file="result/lac.pdf") 


q2 <- ggplot(acid1,aes(time,ace_mean,color=Inoculate))+
  geom_point()+
  geom_line(aes(group=Inoculate))+
  geom_errorbar(aes(ymin=ace_mean-ace_sd,ymax=ace_mean+ace_sd),width = .2)+
  scale_color_npg()+
  scale_y_continuous(limits = c(0,4),n.breaks = 6)+
  theme_bw(base_size = 8)+theme(panel.grid = element_blank(),
                                legend.position = "right",
                                legend.text = element_text(color="black",size=8,face="italic"), # 定义图例文本
                                legend.justification =c(1,0.5) )+
  facet_wrap(Inoculate~.)+
  geom_text(aes(label = sig_ace, y = 3.8),color = "black") +
  xlab(NULL)+ylab("Acetic acid content (mg/g)")+
  theme(axis.title.y = element_text(size = 14),
        strip.text = element_text(size = rel(1.5),face = "italic"),
        strip.background = element_rect(fill = "lightblue", colour ="black",
                                        size = 1)
  )
q2
ggsave(q2,device = "png",file="result/acetic.png") 
ggsave(q2,device = "pdf",file="result/acetic.pdf") 


q3 <- ggplot(acid1,aes(time,eth_mean,color=Inoculate))+
  geom_point()+
  geom_line(aes(group=Inoculate))+
  geom_errorbar(aes(ymin=eth_mean-eth_sd,ymax=eth_mean+eth_sd),width = .2)+
  scale_color_npg()+
  scale_y_continuous(limits = c(0,8),n.breaks = 6)+
  theme_bw(base_size = 8)+theme(panel.grid = element_blank(),
                                legend.position = "right",
                                legend.text = element_text(color="black",size=8,face="italic"), # 定义图例文本
                                legend.justification =c(1,0.5) )+
  facet_wrap(Inoculate~.)+
  geom_text(aes(label = sig_eth, y = 7),color = "black") +
  xlab(NULL)+ylab("Ethanol content (mg/g)")+
  theme(axis.title.y = element_text(size = 14),
        strip.text = element_text(size = rel(1.5),face = "italic"),
        strip.background = element_rect(fill = "lightblue", colour ="black",
                                        size = 1)
  )
q3
ggsave(q3,device = "png",file="result/Ethanol.png") 
ggsave(q3,device = "pdf",file="result/Ethanol.pdf") 


q4 <- ggplot(acid1,aes(time,glu_mean,color=Inoculate))+
  geom_point()+
  geom_line(aes(group=Inoculate))+
  geom_errorbar(aes(ymin=glu_mean-glu_sd,ymax=glu_mean+glu_sd),width = .2)+
  scale_color_npg()+
  scale_y_continuous(limits = c(0,5),n.breaks = 6)+
  theme_bw(base_size = 8)+theme(panel.grid = element_blank(),
                                legend.position = "right",
                                legend.text = element_text(color="black",size=8,face="italic"), # 定义图例文本
                                legend.justification =c(1,0.5) )+
  facet_wrap(Inoculate~.)+
  geom_text(aes(label = sig_glu, y = 4.5),color = "black") +
  xlab(NULL)+ylab("Glucose content (mg/g)")+
  theme(axis.title.y = element_text(size = 14),
        strip.text = element_text(size = rel(1.5),face = "italic"),
        strip.background = element_rect(fill = "lightblue", colour ="black",
                                        size = 1)
  )
q4
ggsave(q4,device = "png",file="result/Glucose.png") 
ggsave(q4,device = "pdf",file="result/Glucose.pdf") 


library(patchwork)
q1 <- q1 + 
  theme(axis.title.y = element_text(size = 16),
        strip.text = element_text(size = rel(1.2),face = "italic"),
        strip.background = element_rect(fill = "lightblue", colour ="black",
                                        size = 1)
  )
q2 <- q2 + 
  theme(axis.title.y = element_text(size = 16),
        strip.text = element_text(size = rel(1.2),face = "italic"),
        strip.background = element_rect(fill = "lightblue", colour ="black",
                                        size = 1)
  )
q3 <- q3 + 
  theme(axis.title.y = element_text(size = 16),
        strip.text = element_text(size = rel(1.2),face = "italic"),
        strip.background = element_rect(fill = "lightblue", colour ="black",
                                        size = 1)
  )
q4 <- q4 + 
  theme(axis.title.y = element_text(size = 16),
        strip.text = element_text(size = rel(1.2),face = "italic"),
        strip.background = element_rect(fill = "lightblue", colour ="black",
                                        size = 1)
  )
a1 <- (q1+q2+q3+q4)+plot_layout(guides = 'collect')
a1
ggsave(a1 ,file="result/yexiang.pdf",width = 24 , height = 14,units = "cm")
ggsave(a1 ,file="result/yexiang.png",width = 24 , height = 14,units = "cm",dpi = 300)


library(cowplot)
q5 <- q1 +
  theme(legend.title=element_text(size=12),legend.text = element_text(size=16))
legend <- get_legend(q5)
mytheme <-  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12))+
  theme(axis.title.y= element_text(size=18))+theme(axis.title.x = element_text(size = 12))+
  theme(legend.title=element_text(size=12),legend.text = element_text(size=12))+theme(legend.position = "none")
p1 <- q1 + mytheme 
p2 <- q2 + mytheme
p3 <- q3 + mytheme 
p4 <- q4 + mytheme


all_1 <- plot_grid(p3, p1 ,p4, p2,  ncol = 1, 
                   align="v", axis = "l",
                   rel_widths = c(1, 1),
                   # labels="AUTO", # auto小写则标签小写
                   # label_size = 18, # 标签大小
                   # label_fontface = "", # 标签字体样式
                   # label_colour = "black", # 标签字体颜色
                   label_x = 0, # 默认为0
                   label_y = 1.0 # 默认为1
) 
#cowp <- plot_grid(all_1, legend, rel_widths = c(3, .6))

ggsave(all_1 ,file="result/yexiang-1_rev.pdf",width = 20 , height = 12,units = "cm")
ggsave(all_1 ,file="result/yexiang-1_rev.png",dpi = 600,width = 10 , height = 15)


library(patchwork)
p1 <- q1 + 
  theme(axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = rel(1.5),face = "italic"),
        strip.background = element_rect(fill = "lightblue", colour ="black",
                                        size = 1)
  )
p2 <- q2 + 
  theme(axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = rel(1.5),face = "italic"),
        strip.background = element_rect(fill = "lightblue", colour ="black",
                                        size = 1)
  )
p3 <- q3 + 
  theme(axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = rel(1.5),face = "italic"),
        strip.background = element_rect(fill = "lightblue", colour ="black",
                                        size = 1)
  )
p4 <- q4 + 
  theme(axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = rel(1.5),face = "italic"),
        strip.background = element_rect(fill = "lightblue", colour ="black",
                                        size = 1)
  )
a1 <- (p1+p2+p3+p4)+plot_layout(guides = 'collect')&
  theme(legend.position='bottom')
a1
ggsave(a1 ,file="result/yx.pdf",width = 24 , height = 14,units = "cm")
ggsave(a1 ,file="result/yx.png",width = 24 , height = 14,units = "cm",dpi = 300)
