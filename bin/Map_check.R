library(ggplot2)

#DEPRECATED
#data <- read.csv("1.2practiques/MCF.Samtools.csv", sep = "\t")
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


TruePositives = Tp_checking("BT_474","edgreen")

names(TruePositives) <- "Fusion"

dataM$Result <- "Novel"
dataM$Result[dataM$Fusion %in% generalM$Gene] <- "Known"
dataM$Result[dataM$Fusion %in% TruePositives$Fusion] <- "TP"


ggplot(data, aes(x= log(pos1), fill=Result)) + geom_boxplot() 
ggplot(data, aes(x= log(pos2), fill=Result)) + geom_boxplot() 


