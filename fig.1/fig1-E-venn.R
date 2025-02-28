#venn-----------------------
library(ggvenn)
library(tidyverse)
library(ggtext)
library(VennDiagram)
co_occ <- read.csv("co_occur.txt" , header=FALSE)
major <-  read.csv("major.txt" , header=FALSE)
ami <-  read.csv("list_aminor.txt" , header=FALSE)
fla <-  read.csv("ist_flaovr.txt" , header=FALSE)


colnames(co_occ) <- "Co_occurring"
colnames(major) <- "Major"
colnames(ami) <- "Amino_acid"
colnames(fla) <- "Flavor"

library(venn)
alist1 <- list(co_occ,major,fla,ami)
alist <- list(Major = major$Major,Co_occurring = co_occ$Co_occurring,
              Amino_acid = ami$Amino_acid,Flavor = fla$Flavor)
alist
venn(alist1)


venn(alist,col = "red",zcolor = "blue")
venn(alist,col = c("red","blue"),zcolor = c("blue","green"))
venn(alist[1:4],col = c("red","blue"),zcolor = c("blue","green"),ellipse = T)
venn(alist,zcolor = rainbow(5),ellipse = T,ilabels =T )


v3 <- venn.diagram(x = alist, filename = NULL, 
                   height = 600, 
                   width = 450,
                   resolution =300, 
                   #imagetype="png", 
                   col="transparent",
                   fill=c("cornflowerblue","green","yellow","darkorchid1"),
                   alpha = 0.50, 
                   cex=1.5, 
                   cat.cex=2
)
cowplot::plot_grid(v3)


