library(tidyverse)
library(vegan)
library(betapart)

group <- read_csv("../data/group_zgj.csv") %>% select(-Repeat)


bac_otu <- read_tsv("../data/contig_taxa_table_species_RC.tsv") %>% column_to_rownames("species") %>% select(-taxonomy)
bac_otu <- t(bac_otu)
bac_otu <- bac_otu %>% decostand(method = "pa")

bac_betapart <- beta.pair(bac_otu,"jaccard")


#all------------------------------------------

bac_betapart_jac <- bac_betapart$beta.jac %>% as.matrix %>% data.frame
bac_betapart_jtu <- bac_betapart$beta.jtu %>% as.matrix %>% data.frame
bac_betapart_jne <- bac_betapart$beta.jne %>% as.matrix %>% data.frame

diag(bac_betapart_jac) <- NA;bac_betapart_jac[lower.tri(bac_betapart_jac)] <- NA
diag(bac_betapart_jtu) <- NA;bac_betapart_jtu[lower.tri(bac_betapart_jtu)] <- NA
diag(bac_betapart_jne) <- NA;bac_betapart_jne[lower.tri(bac_betapart_jne)] <- NA

bac_betapart_jac <- bac_betapart_jac %>% rownames_to_column("sample1") %>% 
  pivot_longer(cols = -1,names_to="sample2",values_to="value") %>% filter(!is.na(value)) %>% add_column(Type="Total")
bac_betapart_jtu <- bac_betapart_jtu %>% rownames_to_column("sample1") %>% 
  pivot_longer(cols = -1,names_to="sample2",values_to="value") %>% filter(!is.na(value)) %>% add_column(Type="Turnover")
bac_betapart_jne <- bac_betapart_jne %>% rownames_to_column("sample1") %>% 
  pivot_longer(cols = -1,names_to="sample2",values_to="value") %>% filter(!is.na(value)) %>% add_column(Type="Nestedness")

bac_betapart_result <- bind_rows(bac_betapart_jac,bac_betapart_jtu,bac_betapart_jne)
bac_betapart_result1 <- left_join(bac_betapart_result,group,by=c("sample1"="Sampleid")) %>% left_join(group,by=c("sample2"="Sampleid"),suffix=c("1","2"))
bac_betapart_result1 <- bac_betapart_result1 %>% select(1:4,SampleTime1,SampleTime2)

bac_betapart_result1$Type <- factor(bac_betapart_result1$Type,levels = c("Total","Turnover","Nestedness"))

bac_betapart_result1$SampleTime1 <- factor(bac_betapart_result1$SampleTime1,levels = unique(group$SampleTime))
bac_betapart_result1$SampleTime2 <- factor(bac_betapart_result1$SampleTime2,levels = unique(group$SampleTime)[length(unique(group$SampleTime)):1])

p1 <- ggplot(bac_betapart_result1,aes(sample1,sample2,fill=value))+
  geom_raster()+scale_fill_viridis_c(limits=c(0,1))+facet_grid(SampleTime2~Type+SampleTime1, scales = "free",switch = "y")+
  theme(strip.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  labs(y="Bacteria SampleTime")


adonis(bac_betapart$beta.jac~SampleTime,group)$aov.tab
adonis(bac_betapart$beta.jtu~SampleTime,group)$aov.tab
adonis(bac_betapart$beta.jne~SampleTime,group)$aov.tab




#CK------------------------------------------------------
group_fer <- group %>% filter(FermentationType=="CK")



bac_otu_fer <- bac_otu[1:15,]
bac_betapart_fer <- beta.pair(bac_otu_fer,"jaccard")

bac_betapart_fer_jac <- bac_betapart_fer$beta.jac %>% as.matrix %>% data.frame
bac_betapart_fer_jtu <- bac_betapart_fer$beta.jtu %>% as.matrix %>% data.frame
bac_betapart_fer_jne <- bac_betapart_fer$beta.jne %>% as.matrix %>% data.frame

diag(bac_betapart_fer_jac) <- NA;bac_betapart_fer_jac[lower.tri(bac_betapart_fer_jac)] <- NA
diag(bac_betapart_fer_jtu) <- NA;bac_betapart_fer_jtu[lower.tri(bac_betapart_fer_jtu)] <- NA
diag(bac_betapart_fer_jne) <- NA;bac_betapart_fer_jne[lower.tri(bac_betapart_fer_jne)] <- NA

bac_betapart_fer_jac <- bac_betapart_fer_jac %>% rownames_to_column("sample1") %>% 
  pivot_longer(cols = -1,names_to="sample2",values_to="value") %>% filter(!is.na(value)) %>% add_column(Type="Total")
bac_betapart_fer_jtu <- bac_betapart_fer_jtu %>% rownames_to_column("sample1") %>% 
  pivot_longer(cols = -1,names_to="sample2",values_to="value") %>% filter(!is.na(value)) %>% add_column(Type="Turnover")
bac_betapart_fer_jne <- bac_betapart_fer_jne %>% rownames_to_column("sample1") %>% 
  pivot_longer(cols = -1,names_to="sample2",values_to="value") %>% filter(!is.na(value)) %>% add_column(Type="Nestedness")

bac_betapart_fer_result <- bind_rows(bac_betapart_fer_jac,bac_betapart_fer_jtu,bac_betapart_fer_jne)
bac_betapart_fer_result <- left_join(bac_betapart_fer_result,group_fer,by=c("sample1"="Sampleid")) %>% left_join(group_fer,by=c("sample2"="Sampleid"),suffix=c("1","2"))
bac_betapart_fer_result <- bac_betapart_fer_result %>% select(1:4,SampleTime1,SampleTime2)

bac_betapart_fer_result$Type <- factor(bac_betapart_fer_result$Type,levels = c("Total","Turnover","Nestedness"))

bac_betapart_fer_result$SampleTime1 <- factor(bac_betapart_fer_result$SampleTime1,levels = unique(group_fer$SampleTime))
bac_betapart_fer_result$SampleTime2 <- factor(bac_betapart_fer_result$SampleTime2,levels = unique(group_fer$SampleTime)[length(unique(group_fer$SampleTime)):1])

p2 <- ggplot(bac_betapart_fer_result,aes(sample1,sample2,fill=value))+
  geom_raster()+scale_fill_viridis_c(limits=c(0,1))+facet_grid(SampleTime2~Type+SampleTime1, scales = "free",switch = "y")+
  theme(strip.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  labs(y="Bacteria SampleTime")


adonis(bac_betapart$beta.jac~SampleTime,group)$aov.tab
adonis(bac_betapart$beta.jtu~SampleTime,group)$aov.tab
adonis(bac_betapart$beta.jne~SampleTime,group)$aov.tab


#Lpscw---------------------------------------------------------
group_fer <- group %>% filter(FermentationType=="Lpscw")



bac_otu_fer <- bac_otu[16:30,]
bac_betapart_fer <- beta.pair(bac_otu_fer,"jaccard")

bac_betapart_fer_jac <- bac_betapart_fer$beta.jac %>% as.matrix %>% data.frame
bac_betapart_fer_jtu <- bac_betapart_fer$beta.jtu %>% as.matrix %>% data.frame
bac_betapart_fer_jne <- bac_betapart_fer$beta.jne %>% as.matrix %>% data.frame

diag(bac_betapart_fer_jac) <- NA;bac_betapart_fer_jac[lower.tri(bac_betapart_fer_jac)] <- NA
diag(bac_betapart_fer_jtu) <- NA;bac_betapart_fer_jtu[lower.tri(bac_betapart_fer_jtu)] <- NA
diag(bac_betapart_fer_jne) <- NA;bac_betapart_fer_jne[lower.tri(bac_betapart_fer_jne)] <- NA

bac_betapart_fer_jac <- bac_betapart_fer_jac %>% rownames_to_column("sample1") %>% 
  pivot_longer(cols = -1,names_to="sample2",values_to="value") %>% filter(!is.na(value)) %>% add_column(Type="Total")
bac_betapart_fer_jtu <- bac_betapart_fer_jtu %>% rownames_to_column("sample1") %>% 
  pivot_longer(cols = -1,names_to="sample2",values_to="value") %>% filter(!is.na(value)) %>% add_column(Type="Turnover")
bac_betapart_fer_jne <- bac_betapart_fer_jne %>% rownames_to_column("sample1") %>% 
  pivot_longer(cols = -1,names_to="sample2",values_to="value") %>% filter(!is.na(value)) %>% add_column(Type="Nestedness")

bac_betapart_fer_result <- bind_rows(bac_betapart_fer_jac,bac_betapart_fer_jtu,bac_betapart_fer_jne)
bac_betapart_fer_result <- left_join(bac_betapart_fer_result,group_fer,by=c("sample1"="Sampleid")) %>% left_join(group_fer,by=c("sample2"="Sampleid"),suffix=c("1","2"))
bac_betapart_fer_result <- bac_betapart_fer_result %>% select(1:4,SampleTime1,SampleTime2)

bac_betapart_fer_result$Type <- factor(bac_betapart_fer_result$Type,levels = c("Total","Turnover","Nestedness"))

bac_betapart_fer_result$SampleTime1 <- factor(bac_betapart_fer_result$SampleTime1,levels = unique(group_fer$SampleTime))
bac_betapart_fer_result$SampleTime2 <- factor(bac_betapart_fer_result$SampleTime2,levels = unique(group_fer$SampleTime)[length(unique(group_fer$SampleTime)):1])

p3 <- ggplot(bac_betapart_fer_result,aes(sample1,sample2,fill=value))+
  geom_raster()+scale_fill_viridis_c(limits=c(0,1))+facet_grid(SampleTime2~Type+SampleTime1, scales = "free",switch = "y")+
  theme(strip.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  labs(y="Bacteria SampleTime")


adonis(bac_betapart$beta.jac~SampleTime,group)$aov.tab
adonis(bac_betapart$beta.jtu~SampleTime,group)$aov.tab
adonis(bac_betapart$beta.jne~SampleTime,group)$aov.tab



#bac----------------------------------------------------------

bac_otu <- read.csv("../data/bac-fungal/bac.csv",row.names = 1) %>% column_to_rownames("taxonomy")
bac_otu <- t(bac_otu)
bac_otu <- bac_otu %>% decostand(method = "pa")

bac_betapart <- beta.pair(bac_otu,"jaccard")


#all------------------------------------------

bac_betapart_jac <- bac_betapart$beta.jac %>% as.matrix %>% data.frame
bac_betapart_jtu <- bac_betapart$beta.jtu %>% as.matrix %>% data.frame
bac_betapart_jne <- bac_betapart$beta.jne %>% as.matrix %>% data.frame

diag(bac_betapart_jac) <- NA;bac_betapart_jac[lower.tri(bac_betapart_jac)] <- NA
diag(bac_betapart_jtu) <- NA;bac_betapart_jtu[lower.tri(bac_betapart_jtu)] <- NA
diag(bac_betapart_jne) <- NA;bac_betapart_jne[lower.tri(bac_betapart_jne)] <- NA

bac_betapart_jac <- bac_betapart_jac %>% rownames_to_column("sample1") %>% 
  pivot_longer(cols = -1,names_to="sample2",values_to="value") %>% filter(!is.na(value)) %>% add_column(Type="Total")
bac_betapart_jtu <- bac_betapart_jtu %>% rownames_to_column("sample1") %>% 
  pivot_longer(cols = -1,names_to="sample2",values_to="value") %>% filter(!is.na(value)) %>% add_column(Type="Turnover")
bac_betapart_jne <- bac_betapart_jne %>% rownames_to_column("sample1") %>% 
  pivot_longer(cols = -1,names_to="sample2",values_to="value") %>% filter(!is.na(value)) %>% add_column(Type="Nestedness")

bac_betapart_result <- bind_rows(bac_betapart_jac,bac_betapart_jtu,bac_betapart_jne)
bac_betapart_result1 <- left_join(bac_betapart_result,group,by=c("sample1"="Sampleid")) %>% left_join(group,by=c("sample2"="Sampleid"),suffix=c("1","2"))
bac_betapart_result1 <- bac_betapart_result1 %>% select(1:4,SampleTime1,SampleTime2)

bac_betapart_result1$Type <- factor(bac_betapart_result1$Type,levels = c("Total","Turnover","Nestedness"))

bac_betapart_result1$SampleTime1 <- factor(bac_betapart_result1$SampleTime1,levels = unique(group$SampleTime))
bac_betapart_result1$SampleTime2 <- factor(bac_betapart_result1$SampleTime2,levels = unique(group$SampleTime)[length(unique(group$SampleTime)):1])

p4 <- ggplot(bac_betapart_result1,aes(sample1,sample2,fill=value))+
  geom_raster()+scale_fill_viridis_c(limits=c(0,1))+facet_grid(SampleTime2~Type+SampleTime1, scales = "free",switch = "y")+
  theme(strip.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  labs(y="Bacteria SampleTime")


adonis(bac_betapart$beta.jac~SampleTime,group)$aov.tab
adonis(bac_betapart$beta.jtu~SampleTime,group)$aov.tab
adonis(bac_betapart$beta.jne~SampleTime,group)$aov.tab




#CK------------------------------------------------------
group_fer <- group %>% filter(FermentationType=="CK")



bac_otu_fer <- bac_otu[1:15,]
bac_betapart_fer <- beta.pair(bac_otu_fer,"jaccard")

bac_betapart_fer_jac <- bac_betapart_fer$beta.jac %>% as.matrix %>% data.frame
bac_betapart_fer_jtu <- bac_betapart_fer$beta.jtu %>% as.matrix %>% data.frame
bac_betapart_fer_jne <- bac_betapart_fer$beta.jne %>% as.matrix %>% data.frame

diag(bac_betapart_fer_jac) <- NA;bac_betapart_fer_jac[lower.tri(bac_betapart_fer_jac)] <- NA
diag(bac_betapart_fer_jtu) <- NA;bac_betapart_fer_jtu[lower.tri(bac_betapart_fer_jtu)] <- NA
diag(bac_betapart_fer_jne) <- NA;bac_betapart_fer_jne[lower.tri(bac_betapart_fer_jne)] <- NA

bac_betapart_fer_jac <- bac_betapart_fer_jac %>% rownames_to_column("sample1") %>% 
  pivot_longer(cols = -1,names_to="sample2",values_to="value") %>% filter(!is.na(value)) %>% add_column(Type="Total")
bac_betapart_fer_jtu <- bac_betapart_fer_jtu %>% rownames_to_column("sample1") %>% 
  pivot_longer(cols = -1,names_to="sample2",values_to="value") %>% filter(!is.na(value)) %>% add_column(Type="Turnover")
bac_betapart_fer_jne <- bac_betapart_fer_jne %>% rownames_to_column("sample1") %>% 
  pivot_longer(cols = -1,names_to="sample2",values_to="value") %>% filter(!is.na(value)) %>% add_column(Type="Nestedness")

bac_betapart_fer_result <- bind_rows(bac_betapart_fer_jac,bac_betapart_fer_jtu,bac_betapart_fer_jne)
bac_betapart_fer_result <- left_join(bac_betapart_fer_result,group_fer,by=c("sample1"="Sampleid")) %>% left_join(group_fer,by=c("sample2"="Sampleid"),suffix=c("1","2"))
bac_betapart_fer_result <- bac_betapart_fer_result %>% select(1:4,SampleTime1,SampleTime2)

bac_betapart_fer_result$Type <- factor(bac_betapart_fer_result$Type,levels = c("Total","Turnover","Nestedness"))

bac_betapart_fer_result$SampleTime1 <- factor(bac_betapart_fer_result$SampleTime1,levels = unique(group_fer$SampleTime))
bac_betapart_fer_result$SampleTime2 <- factor(bac_betapart_fer_result$SampleTime2,levels = unique(group_fer$SampleTime)[length(unique(group_fer$SampleTime)):1])

p5 <- ggplot(bac_betapart_fer_result,aes(sample1,sample2,fill=value))+
  geom_raster()+scale_fill_viridis_c(limits=c(0,1))+facet_grid(SampleTime2~Type+SampleTime1, scales = "free",switch = "y")+
  theme(strip.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  labs(y="Bacteria SampleTime")


adonis(bac_betapart$beta.jac~SampleTime,group)$aov.tab
adonis(bac_betapart$beta.jtu~SampleTime,group)$aov.tab
adonis(bac_betapart$beta.jne~SampleTime,group)$aov.tab


#Lpscw---------------------------------------------------------
group_fer <- group %>% filter(FermentationType=="Lpscw")



bac_otu_fer <- bac_otu[16:30,]
bac_betapart_fer <- beta.pair(bac_otu_fer,"jaccard")

bac_betapart_fer_jac <- bac_betapart_fer$beta.jac %>% as.matrix %>% data.frame
bac_betapart_fer_jtu <- bac_betapart_fer$beta.jtu %>% as.matrix %>% data.frame
bac_betapart_fer_jne <- bac_betapart_fer$beta.jne %>% as.matrix %>% data.frame

diag(bac_betapart_fer_jac) <- NA;bac_betapart_fer_jac[lower.tri(bac_betapart_fer_jac)] <- NA
diag(bac_betapart_fer_jtu) <- NA;bac_betapart_fer_jtu[lower.tri(bac_betapart_fer_jtu)] <- NA
diag(bac_betapart_fer_jne) <- NA;bac_betapart_fer_jne[lower.tri(bac_betapart_fer_jne)] <- NA

bac_betapart_fer_jac <- bac_betapart_fer_jac %>% rownames_to_column("sample1") %>% 
  pivot_longer(cols = -1,names_to="sample2",values_to="value") %>% filter(!is.na(value)) %>% add_column(Type="Total")
bac_betapart_fer_jtu <- bac_betapart_fer_jtu %>% rownames_to_column("sample1") %>% 
  pivot_longer(cols = -1,names_to="sample2",values_to="value") %>% filter(!is.na(value)) %>% add_column(Type="Turnover")
bac_betapart_fer_jne <- bac_betapart_fer_jne %>% rownames_to_column("sample1") %>% 
  pivot_longer(cols = -1,names_to="sample2",values_to="value") %>% filter(!is.na(value)) %>% add_column(Type="Nestedness")

bac_betapart_fer_result <- bind_rows(bac_betapart_fer_jac,bac_betapart_fer_jtu,bac_betapart_fer_jne)
bac_betapart_fer_result <- left_join(bac_betapart_fer_result,group_fer,by=c("sample1"="Sampleid")) %>% left_join(group_fer,by=c("sample2"="Sampleid"),suffix=c("1","2"))
bac_betapart_fer_result <- bac_betapart_fer_result %>% select(1:4,SampleTime1,SampleTime2)

bac_betapart_fer_result$Type <- factor(bac_betapart_fer_result$Type,levels = c("Total","Turnover","Nestedness"))

bac_betapart_fer_result$SampleTime1 <- factor(bac_betapart_fer_result$SampleTime1,levels = unique(group_fer$SampleTime))
bac_betapart_fer_result$SampleTime2 <- factor(bac_betapart_fer_result$SampleTime2,levels = unique(group_fer$SampleTime)[length(unique(group_fer$SampleTime)):1])

p6 <- ggplot(bac_betapart_fer_result,aes(sample1,sample2,fill=value))+
  geom_raster()+scale_fill_viridis_c(limits=c(0,1))+facet_grid(SampleTime2~Type+SampleTime1, scales = "free",switch = "y")+
  theme(strip.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  labs(y="Bacteria SampleTime")


adonis(bac_betapart$beta.jac~SampleTime,group)$aov.tab
adonis(bac_betapart$beta.jtu~SampleTime,group)$aov.tab
adonis(bac_betapart$beta.jne~SampleTime,group)$aov.tab


#Enka----------------------------------------------------------

bac_otu <- read.csv("../data/bac-fungal/Euka.csv",row.names = 1) %>% column_to_rownames("taxonomy")
bac_otu <- t(bac_otu)
bac_otu <- bac_otu %>% decostand(method = "pa")

bac_betapart <- beta.pair(bac_otu,"jaccard")


#all------------------------------------------

bac_betapart_jac <- bac_betapart$beta.jac %>% as.matrix %>% data.frame
bac_betapart_jtu <- bac_betapart$beta.jtu %>% as.matrix %>% data.frame
bac_betapart_jne <- bac_betapart$beta.jne %>% as.matrix %>% data.frame

diag(bac_betapart_jac) <- NA;bac_betapart_jac[lower.tri(bac_betapart_jac)] <- NA
diag(bac_betapart_jtu) <- NA;bac_betapart_jtu[lower.tri(bac_betapart_jtu)] <- NA
diag(bac_betapart_jne) <- NA;bac_betapart_jne[lower.tri(bac_betapart_jne)] <- NA

bac_betapart_jac <- bac_betapart_jac %>% rownames_to_column("sample1") %>% 
  pivot_longer(cols = -1,names_to="sample2",values_to="value") %>% filter(!is.na(value)) %>% add_column(Type="Total")
bac_betapart_jtu <- bac_betapart_jtu %>% rownames_to_column("sample1") %>% 
  pivot_longer(cols = -1,names_to="sample2",values_to="value") %>% filter(!is.na(value)) %>% add_column(Type="Turnover")
bac_betapart_jne <- bac_betapart_jne %>% rownames_to_column("sample1") %>% 
  pivot_longer(cols = -1,names_to="sample2",values_to="value") %>% filter(!is.na(value)) %>% add_column(Type="Nestedness")

bac_betapart_result <- bind_rows(bac_betapart_jac,bac_betapart_jtu,bac_betapart_jne)
bac_betapart_result1 <- left_join(bac_betapart_result,group,by=c("sample1"="Sampleid")) %>% left_join(group,by=c("sample2"="Sampleid"),suffix=c("1","2"))
bac_betapart_result1 <- bac_betapart_result1 %>% select(1:4,SampleTime1,SampleTime2)

bac_betapart_result1$Type <- factor(bac_betapart_result1$Type,levels = c("Total","Turnover","Nestedness"))

bac_betapart_result1$SampleTime1 <- factor(bac_betapart_result1$SampleTime1,levels = unique(group$SampleTime))
bac_betapart_result1$SampleTime2 <- factor(bac_betapart_result1$SampleTime2,levels = unique(group$SampleTime)[length(unique(group$SampleTime)):1])

p7 <- ggplot(bac_betapart_result1,aes(sample1,sample2,fill=value))+
  geom_raster()+scale_fill_viridis_c(limits=c(0,1))+facet_grid(SampleTime2~Type+SampleTime1, scales = "free",switch = "y")+
  theme(strip.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  labs(y="Bacteria SampleTime")


adonis(bac_betapart$beta.jac~SampleTime,group)$aov.tab
adonis(bac_betapart$beta.jtu~SampleTime,group)$aov.tab
adonis(bac_betapart$beta.jne~SampleTime,group)$aov.tab




#CK------------------------------------------------------
group_fer <- group %>% filter(FermentationType=="CK")



bac_otu_fer <- bac_otu[1:15,]
bac_betapart_fer <- beta.pair(bac_otu_fer,"jaccard")

bac_betapart_fer_jac <- bac_betapart_fer$beta.jac %>% as.matrix %>% data.frame
bac_betapart_fer_jtu <- bac_betapart_fer$beta.jtu %>% as.matrix %>% data.frame
bac_betapart_fer_jne <- bac_betapart_fer$beta.jne %>% as.matrix %>% data.frame

diag(bac_betapart_fer_jac) <- NA;bac_betapart_fer_jac[lower.tri(bac_betapart_fer_jac)] <- NA
diag(bac_betapart_fer_jtu) <- NA;bac_betapart_fer_jtu[lower.tri(bac_betapart_fer_jtu)] <- NA
diag(bac_betapart_fer_jne) <- NA;bac_betapart_fer_jne[lower.tri(bac_betapart_fer_jne)] <- NA

bac_betapart_fer_jac <- bac_betapart_fer_jac %>% rownames_to_column("sample1") %>% 
  pivot_longer(cols = -1,names_to="sample2",values_to="value") %>% filter(!is.na(value)) %>% add_column(Type="Total")
bac_betapart_fer_jtu <- bac_betapart_fer_jtu %>% rownames_to_column("sample1") %>% 
  pivot_longer(cols = -1,names_to="sample2",values_to="value") %>% filter(!is.na(value)) %>% add_column(Type="Turnover")
bac_betapart_fer_jne <- bac_betapart_fer_jne %>% rownames_to_column("sample1") %>% 
  pivot_longer(cols = -1,names_to="sample2",values_to="value") %>% filter(!is.na(value)) %>% add_column(Type="Nestedness")

bac_betapart_fer_result <- bind_rows(bac_betapart_fer_jac,bac_betapart_fer_jtu,bac_betapart_fer_jne)
bac_betapart_fer_result <- left_join(bac_betapart_fer_result,group_fer,by=c("sample1"="Sampleid")) %>% left_join(group_fer,by=c("sample2"="Sampleid"),suffix=c("1","2"))
bac_betapart_fer_result <- bac_betapart_fer_result %>% select(1:4,SampleTime1,SampleTime2)

bac_betapart_fer_result$Type <- factor(bac_betapart_fer_result$Type,levels = c("Total","Turnover","Nestedness"))

bac_betapart_fer_result$SampleTime1 <- factor(bac_betapart_fer_result$SampleTime1,levels = unique(group_fer$SampleTime))
bac_betapart_fer_result$SampleTime2 <- factor(bac_betapart_fer_result$SampleTime2,levels = unique(group_fer$SampleTime)[length(unique(group_fer$SampleTime)):1])

p8 <- ggplot(bac_betapart_fer_result,aes(sample1,sample2,fill=value))+
  geom_raster()+scale_fill_viridis_c(limits=c(0,1))+facet_grid(SampleTime2~Type+SampleTime1, scales = "free",switch = "y")+
  theme(strip.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  labs(y="Bacteria SampleTime")


adonis(bac_betapart$beta.jac~SampleTime,group)$aov.tab
adonis(bac_betapart$beta.jtu~SampleTime,group)$aov.tab
adonis(bac_betapart$beta.jne~SampleTime,group)$aov.tab


#Lpscw---------------------------------------------------------
group_fer <- group %>% filter(FermentationType=="Lpscw")



bac_otu_fer <- bac_otu[16:30,]
bac_betapart_fer <- beta.pair(bac_otu_fer,"jaccard")

bac_betapart_fer_jac <- bac_betapart_fer$beta.jac %>% as.matrix %>% data.frame
bac_betapart_fer_jtu <- bac_betapart_fer$beta.jtu %>% as.matrix %>% data.frame
bac_betapart_fer_jne <- bac_betapart_fer$beta.jne %>% as.matrix %>% data.frame

diag(bac_betapart_fer_jac) <- NA;bac_betapart_fer_jac[lower.tri(bac_betapart_fer_jac)] <- NA
diag(bac_betapart_fer_jtu) <- NA;bac_betapart_fer_jtu[lower.tri(bac_betapart_fer_jtu)] <- NA
diag(bac_betapart_fer_jne) <- NA;bac_betapart_fer_jne[lower.tri(bac_betapart_fer_jne)] <- NA

bac_betapart_fer_jac <- bac_betapart_fer_jac %>% rownames_to_column("sample1") %>% 
  pivot_longer(cols = -1,names_to="sample2",values_to="value") %>% filter(!is.na(value)) %>% add_column(Type="Total")
bac_betapart_fer_jtu <- bac_betapart_fer_jtu %>% rownames_to_column("sample1") %>% 
  pivot_longer(cols = -1,names_to="sample2",values_to="value") %>% filter(!is.na(value)) %>% add_column(Type="Turnover")
bac_betapart_fer_jne <- bac_betapart_fer_jne %>% rownames_to_column("sample1") %>% 
  pivot_longer(cols = -1,names_to="sample2",values_to="value") %>% filter(!is.na(value)) %>% add_column(Type="Nestedness")

bac_betapart_fer_result <- bind_rows(bac_betapart_fer_jac,bac_betapart_fer_jtu,bac_betapart_fer_jne)
bac_betapart_fer_result <- left_join(bac_betapart_fer_result,group_fer,by=c("sample1"="Sampleid")) %>% left_join(group_fer,by=c("sample2"="Sampleid"),suffix=c("1","2"))
bac_betapart_fer_result <- bac_betapart_fer_result %>% select(1:4,SampleTime1,SampleTime2)

bac_betapart_fer_result$Type <- factor(bac_betapart_fer_result$Type,levels = c("Total","Turnover","Nestedness"))

bac_betapart_fer_result$SampleTime1 <- factor(bac_betapart_fer_result$SampleTime1,levels = unique(group_fer$SampleTime))
bac_betapart_fer_result$SampleTime2 <- factor(bac_betapart_fer_result$SampleTime2,levels = unique(group_fer$SampleTime)[length(unique(group_fer$SampleTime)):1])

p9 <- ggplot(bac_betapart_fer_result,aes(sample1,sample2,fill=value))+
  geom_raster()+scale_fill_viridis_c(limits=c(0,1))+facet_grid(SampleTime2~Type+SampleTime1, scales = "free",switch = "y")+
  theme(strip.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  labs(y="Bacteria SampleTime")+ theme(legend.position="none")


adonis(bac_betapart$beta.jac~SampleTime,group)$aov.tab
adonis(bac_betapart$beta.jtu~SampleTime,group)$aov.tab
adonis(bac_betapart$beta.jne~SampleTime,group)$aov.tab
