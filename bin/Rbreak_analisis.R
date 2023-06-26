library(tidyverse)
library(dplyr)
library(ggplot2)
bks <- read.csv("1.2practiques/BREAK/REAL_SCORES_NODUPS_RSHIFT.csv")
cbks <-bks
bks$Fusion <- as.character(bks$Fusion)
bks$Result <- "Novel"
bks$Result[bks$Fusion %in% generalM$Gene] <- "Known"
bks$Result[bks$Fusion %in% TruePositives$Fusion] <- "TP"
intersect(bks$Fusion,TruePositives$Fusion)
bks
bks <- bks[bks$Gene1Normal != "noR",]
bks$Gene1Normal <- as.numeric(bks$Gene1Normal)
bks$Gene1Fused <- as.numeric(bks$Gene1Fused)
bks$Gene2Normal <- as.numeric(bks$Gene2Normal)
bks$Gene2Fused <- as.numeric(bks$Gene2Fused)
bks$Fused1 <- as.numeric(bks$Fused1)
bks$Fused2 <- as.numeric(bks$Fused2)
bks$Anchor1 <- as.numeric(bks$Anchor1)
bks$Anchor2 <- as.numeric(bks$Anchor2)
bks$Gene1total <- as.numeric(bks$Gene1total)
bks$Gene2total <- as.numeric(bks$Gene2total)
bks %>% dplyr::select(!c(Fusion,Anchor2,Anchor1,Fused2,Fused1,Gene2Fused,Gene1Fused,Gene2Normal,Gene1Normal)) %>%
  pivot_longer(!Result) %>%
  ggplot(aes(x=name,y=value,fill=name))+
  geom_boxplot()+
  facet_wrap(.~Result,scale='free')
  
ggplot(bks,aes(x=Gene1Fused, y=Gene2Fused))+ geom_point()
mean(bks$Gene1total)
mean(bks$Gene2total)
mean(bks$Fused1)
mean(bks$Fused2)
