library(stringr)
library(ggvenn)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
#READ STANDARIZED DATA
dataM <- read.csv("1.2practiques/Standarized_calls.csv")

#READ MITLEMAN DICTIONARIES
MorphDict<- read.csv("1.2practiques/mitelman_db/MITDICTMORPH.TXT", sep = ",",header = FALSE)
TopDict <- read.csv("1.2practiques/mitelman_db/MITDICTTOP.TXT", sep = "," ,header = FALSE)
generalM <- read.csv("1.2practiques/mitelman_db/MCGENE.TXT.DATA", sep = "\t")
Annot <- read.csv("1.2practiques/mitelman_db/MBCA.TXT.DATA", sep = "\t")

dataM$mitleman <- "N"
dataM$mitleman[dataM$Fusion %in% generalM$Gene] <- "Y"

oncogeneicity <- read.csv("1.2practiques/ongene_human.txt", sep = "\t")



#FUNCTIOOONS

translate_morph <- function(morphology){
  if (morphology =="All"){
    # sele <- MorphDict[!grepl("%",MorphDict$V1, fixed = TRUE),]
    return("All")
  }
  sele <- MorphDict[MorphDict$V2 %in% morphology,]
  if (grepl("%", sele$V1, fixed = TRUE)){
      search <- substr(sele$V1, 1,2)
      sele <- MorphDict[substr(MorphDict$V1,1,2) == search,]
      quit <- paste0(search,"%")
      sele <- sele[sele$V1 != quit,]
  }
  return(sele$V1)
}

translate_top <- function(topology){
  if (topology =="All"){
    # sele <- TopDict[!grepl("%",TopDict$V1, fixed = TRUE),]
    return("All")
  }
  sele <- TopDict[TopDict$V2 %in% topology,]
  if (grepl("%", sele$V1, fixed = TRUE)){
    search <- substr(sele$V1, 1,2)
    sele <- TopDict[substr(TopDict$V1,1,2) == search,]
    quit <- paste0(search,"%")
    sele <- sele[sele$V1 != quit,]
  }
  return(sele$V1)
}

unite_choices <- function(choices, f){
  chosen <- c()
  for (i in choices){
    if(f == "T"){
      chosen <- append(chosen, translate_top(i))
    }
    else{
      chosen <- append(chosen, translate_morph(i))
    }
  }
  print(unique(chosen))
  return(unique(chosen))
}

check_mitelman <- function(fusion_str_c, morphologies , topologies){
  topologies <- unite_choices(topologies,"T")
  morphologies <- unite_choices(morphologies,"M")
  print("Computing topology table")
  Top<- check_topology(topologies)
  print("Computing morph table")
  Morph <- check_morphology(morphologies)
  print("Computing union table")
  Union <- check_union(topologies,morphologies)
 
  xlist <- list()
  c = 1
  for (fusion_str in fusion_str_c){
    m = 0
    t = 0 
    u = 0
    if (fusion_str %in% generalM$Gene){
      refe <- generalM$RefNo[generalM$Gene == fusion_str]
      for (ref in refe){
        
        if (nrow(Morph[Morph$RefNo == ref,])>0){
          m = 1
        }
        if (nrow(Top[Top$RefNo == ref,])>0){
          t = 1
        }
        if (nrow(Union[Union$RefNo == ref,])>0){
          u = 1
        }
      }
      xlist[[c]] <- c(1,m,t,u)
    }
    else{
      xlist[[c]] <- c(0,0,0,0)
    }
    c = c +1
  }
  return(xlist)
}

check_union <- function(topologies, morphologies){
  NotNA <- Annot[!is.na(Annot$Morph),]
  NotNA <- NotNA[!is.na(NotNA$Topo),]
  if (topologies == morphologies){
    return(NotNA)
  }
  tot <- data_frame()

  for (i in morphologies){
    sele1 <- NotNA[NotNA$Morph == as.integer(i),]
    for (j in topologies){
      sele <- sele1[sele1$Topo == as.integer(j),]
      tot <- rbind(tot,sele)
      }
      
  }
  return(tot)
}
check_morphology <- function(morphologies){
  NotNA <- Annot[!is.na(Annot$Morph),]
  if (morphologies == "All"){
    return(NotNA)
  }
  tot <- data_frame()
  for (i in morphologies){
      
      sele <- NotNA[NotNA$Morph == as.integer(i),]
      tot <- rbind(tot,sele)
  }
  return(tot)
}

check_topology <- function(topologies){
  tot <- data_frame()
  NotNA <- Annot[!is.na(Annot$Topo),]
  if(topologies == "All"){
    return(NotNA)
  }
  for (i in topologies){
    sele <- NotNA[NotNA$Topo == as.integer(i) ,]
    tot <- rbind(tot, sele)
  }
  return(tot)
}

################################################
#!!! Add to standarization the kind of translocation dup etc check_form <- function()

Tp_checking <- function(sample_name, exp_name){
  Combinations <- c()
  experiment <- read.csv(file = paste0("1.2practiques/",paste0(exp_name,"/fusions.txt")), sep = "\t")
  experiment_sel <- experiment[experiment$Sample == sample_name,]
  for (num in 1:length(experiment_sel$Gene1)){
    A <- experiment_sel$Gene1[num]
    A <- unlist(strsplit(A, ','))
    B <- experiment_sel$Gene2[num]
    B <- unlist(strsplit(B, ','))
    for(i in A){
      for(j in B){
        Combinations<-append(Combinations,paste(i,j,sep="::"))
      }
    }
  }
  return(as.data.frame(Combinations))
}


TruePositives = Tp_checking(sample,"edgreen")

names(TruePositives) <- "Fusion"

dataM$Result <- "Novel"
dataM$Result[dataM$Fusion %in% generalM$Gene] <- "Known"
dataM$Result[dataM$Fusion %in% TruePositives$Fusion] <- "TP"

check_onco <- function(Gene1_c, Gene2_c){
  oncogene <- c()
  for( i in 1:length(Gene1_c)){
    t = 0
    if (sum(grepl(Gene1_c[i],oncogeneicity$OncogeneName)) || sum(grepl(Gene1_c[i],oncogeneicity$Alias))){
      t = t + 1
    }
    if (sum(grepl(Gene2_c[i],oncogeneicity$OncogeneName)) || sum(grepl(Gene2_c[i],oncogeneicity$Alias))){
      t = t + 1
    }
    oncogene <- append(oncogene,t)
  }
  return(oncogene)
}






































df.long <- dataM %>% 
  pivot_longer(spanning.pairs:spanning.reads, names_to = 'Statistic', values_to = 'value') %>% 
  group_by('Statistic')

df.long$o <- df.long$value == 0
tiff("test.tiff", units="in", width=5, height=5, res=300)
ggplot(df.long, aes(x= Result ,value, fill = Statistic)) + geom_boxplot() + coord_flip()  + theme_minimal()  +labs(fill="Statistic", title= "Spanning Reads and Pairs across Fusion Call types") +theme(legend.position = "bottom", plot.title = element_text(family = "Luminari", face = "bold", size = (12)))+ ylab("Amount") + xlab("Fusion Call Type")+ guides(fill = guide_legend(reverse=TRUE)) + scale_fill_discrete(name = "Statistic", labels = c("Spanning Pairs", "Spanning Reads")) + geom_hline(yintercept = 2)

dev.off()

ggplot(dataM, aes(x= spanning.reads, fill=Result)) + geom_boxplot() 

selected <- generalM[grepl("::",generalM$Gene),]
selected_diss = Annot[Annot$RefNo %in% selected$RefNo,]
uniq<- unique(selected_diss$Morph)
morp_sel <- MorphDict[MorphDict$V1 %in% uniq,]
names(morp_sel) <- c("Morph", "MorphN")
morp_sel$RefNo <- as.numeric(morp_sel$Morph)
loco <- base::merge(morp_sel,selected_diss,all.y = TRUE, by="Morph")
tabla<-table(loco$MorphN)
data <- as.data.frame(tabla)
data <- data[order(data$Freq, decreasing = TRUE),]
dataP <- data[1:20,]
calc = sum(data$Freq)-sum(dataP$Freq)
dataP$Var1 <- as.character(dataP$Var1)
dataP[nrow(dataP)+1,] <- list("others", calc )

dataP$Var1
loco$PlotM <- loco$MorphN %in% dataP$Var1
length(loco$MorphN %in% dataP$Var1)
loco$MorphN[loco$PlotM == FALSE] <- "Others"
table(loco$MorphN)

mycolors <- colorRampPalette(brewer.pal(9, "PuBu"))(30)
pt<- ggplot(loco, aes(y = MorphN, fill =MorphN)) + geom_bar() + theme(legend.position = "None", axis.title.y = element_blank(), plot.title = element_text(family = "Luminari", face = "bold", size = (12))) + scale_fill_manual(values = mycolors[9:30])
pt
#### TOPOLOGY

selected_diss <- selected_diss[!is.na(selected_diss$Topo),]
uniq<- unique(selected_diss$Topo)
TopDict$V1 <- as.numeric(TopDict$V1)
top_sel <- TopDict[TopDict$V1 %in% uniq,]

names(top_sel) <- c("Topo", "TopoN")
top_sel$Topo <- as.numeric(top_sel$Topo)
loco <- base::merge(top_sel,selected_diss,all.y = TRUE, by="Topo")
loco$TopoN[loco$Topo == 0] <- "All" 
loco <- loco[loco$TopoN != "Unknown site",]
tabla<-table(loco$TopoN)
data <- as.data.frame(tabla)
data <- data[order(data$Freq, decreasing = TRUE),]
dataP <- data[1:14,]
calc = sum(dataP$Freq)-sum(dataP$Freq)
dataP$Var1 <- as.character(dataP$Var1)
dataP[nrow(dataP)+1,] <- list("Others", calc )

dataP$Var1
length(loco$PlotM)
loco$PlotM <- loco$TopoN %in% dataP$Var1
length(loco$TopoN %in% dataP$Var1)
loco$TopoN[loco$PlotM == FALSE] <- "Others"




pa <- ggplot(loco, aes(x = TopoN, fill = TopoN)) + geom_bar() + labs(x = "Tissue", y= "count" )+ theme(legend.position = "None") + scale_fill_manual(values = mycolors[9:30]) + 
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) 
pa
lko <-ggarrange(pt,pa)
lko
annotate_figure(lko, top = text_grob("Number of anotated fusions by disease and tissue in the Mitelman Database", 
                                      color = "black", face = "bold", size = 17))
##########################################

check_onco(dataM$Gene1,dataM$Gene2)

