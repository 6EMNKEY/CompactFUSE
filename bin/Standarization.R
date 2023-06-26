library(stringr)
library(ggvenn)
library(dplyr)
library(ggplot2)
library("org.Mm.eg.db")
library("org.Hs.eg.db")
library("biomaRt")
??org.Mm.eg.db
#Remember to add option for junk prio??

load_fusioncatcher <- function(exp_name,samp_name){
  #fusion catcher
  fusion_catcher = read.csv(sprintf("1.2practiques/%s/FusionCatcher/%s/final-list_candidate-fusion-genes.txt", exp_name, samp_name),sep = "\t", header = TRUE)
  names(fusion_catcher)[1] = "Gene1"
  names(fusion_catcher)[2] = "Gene2"
  fusion_catcher["Fusion"] = paste(sep="::", fusion_catcher$Gene1, fusion_catcher$Gene2)
  return(fusion_catcher)
}
load_jaffa <- function(exp_name,samp_name){
  #JAFFA
  jaffa = read.csv(sprintf("1.2practiques/%s/JAFFA/%s/jaffa_results.csv", exp_name, samp_name))
  jaffa[c("Gene1", "Gene2")] = str_split_fixed(jaffa$fusion.genes, ':', n= 2)
  jaffa = jaffa[names(jaffa)!= "fusion.genes"]
  jaffa["Fusion"] = paste(sep = "::", jaffa$Gene1, jaffa$Gene2)
  return(jaffa)
}
load_arriba <- function(exp_name,samp_name){
  #ARRIBA
  arriba_tsv = read.csv(sprintf("1.2practiques/%s/Arriba/%s/fusions.tsv", exp_name, samp_name), sep = "\t")
  names(arriba_tsv)[names(arriba_tsv)== "X.gene1"] <- "Gene1"
  names(arriba_tsv)[names(arriba_tsv)== "gene2"] <- "Gene2"
  arriba_tsv["Fusion"] = paste(sep = "::", arriba_tsv$Gene1, arriba_tsv$Gene2)
  return(arriba_tsv)
}
  
load_starFusion <- function(exp_name,samp_name){
  #STAR-FUSION
  star_fusion = read.csv(sprintf("1.2practiques/%s/StarFusion/%s/STAR-Fusion_outdir/star-fusion.fusion_predictions.abridged.tsv", exp_name, samp_name), sep = "\t")
  star_fusion[c("Gene1", "Gene2")] = str_split_fixed(star_fusion$X.FusionName, "--", n = 2)
  star_fusion = star_fusion[names(star_fusion) != "X.FusionName"]
  star_fusion["Fusion"] = paste(sep = "::", star_fusion$Gene1, star_fusion$Gene2)
  return(star_fusion)
}
sample <- "BT_474"
sample2 <-  "BT_474"
fusion_catcher <- ""
jaffa <- ""
arriba_tsv <- ""
star_fusion <- ""
fusion_catcher <- load_fusioncatcher("edgreen", sample2)
jaffa <- load_jaffa("edgreen", sample)
arriba_tsv <- load_arriba("edgreen", sample)
star_fusion <- load_starFusion("edgreen", sample)


standarise <- function(){
  #STANDARIZED
  #JAFFA has no gene id's (have to be added)
  jaffa$chrom1<- substr(jaffa$chrom1, 4,5)
  jaffa$chrom2<- substr(jaffa$chrom2, 4,5)
  K1 <- dplyr::select(jaffa,Fusion,Gene1,Gene2,base1,chrom1,strand1,base2,chrom2,strand2,spanning.pairs,spanning.reads)
  K1$GeneId1 <- "."
  K1$GeneId2 <- "."
  K1$caller <- "JAFFA"

  #Fusion_catcher

  if (length(fusion_catcher$Gene1) >= 1){
    colnames(fusion_catcher)[colnames(fusion_catcher)=="Spanning_pairs"]<- "spanning.pairs"
    colnames(fusion_catcher)[colnames(fusion_catcher)=="Spanning_unique_reads"]<- "spanning.reads"
    fusion_catcher[c("chrom1", "base1","strand1")] = str_split_fixed(fusion_catcher$Fusion_point_for_gene_2.3end_fusion_partner., ':', n= 3)
    fusion_catcher[c("chrom2", "base2","strand2")] = str_split_fixed(fusion_catcher$Fusion_point_for_gene_1.5end_fusion_partner., ':', n= 3)
    K2 <- dplyr::select(fusion_catcher,Fusion,Gene1,Gene2,base1,chrom1,strand1,base2,chrom2,strand2,spanning.pairs,spanning.reads)
    K2$GeneId1 <- fusion_catcher$Gene_2_id.3end_fusion_partner.
    K2$GeneId2 <- fusion_catcher$Gene_1_id.5end_fusion_partner.
    K2$caller <- "FusionCatcher"
  }
  else{
    #K2 <- as.data.frame(o="A",a ="b", b = "c", d ="",row.names = names(K1))
  }
  
  # #Arriba
  colnames(arriba_tsv)[colnames(arriba_tsv)=="discordant_mates"]<- "spanning.reads"
  arriba_tsv$spanning.pairs <-arriba_tsv$split_reads1 +arriba_tsv$split_reads2
  arriba_tsv[c("chrom1", "base1")] = str_split_fixed(arriba_tsv$breakpoint1, ':', n= 2)
  arriba_tsv[c("chrom2", "base2")] = str_split_fixed(arriba_tsv$breakpoint2, ':', n= 2)
  colnames(arriba_tsv)[colnames(arriba_tsv)=="strand1.gene.fusion."]<- "strand1"
  colnames(arriba_tsv)[colnames(arriba_tsv)=="strand2.gene.fusion."]<- "strand2"
  K3 <- dplyr::select(arriba_tsv,Fusion,Gene1,Gene2,base1,chrom1,strand1,base2,chrom2,strand2,spanning.pairs,spanning.reads)
  K3$GeneId1 <- arriba_tsv$gene_id1
  K3$GeneId2 <- arriba_tsv$gene_id2
  K3$caller <- "Arriba"

  #StarFusion

  colnames(star_fusion)[colnames(star_fusion)=="JunctionReadCount"]<- "spanning.reads"
  colnames(star_fusion)[colnames(star_fusion)=="SpanningFragCount"]<- "spanning.pairs"
  star_fusion[c("chrom1", "base1","strand1")] = str_split_fixed(star_fusion$LeftBreakpoint, ':', n= 3)
  star_fusion[c("chrom2", "base2","strand2")] = str_split_fixed(star_fusion$RightBreakpoint, ':', n= 3)
  star_fusion[c("t", "GeneId1")] = str_split_fixed(star_fusion$LeftGene, "\\^", n= 2)
  star_fusion[c("t", "GeneId2")] = str_split_fixed(star_fusion$RightGene, "\\^", n= 2)
  star_fusion$chrom1<- substr(star_fusion$chrom1, 4,5)
  star_fusion$chrom2<- substr(star_fusion$chrom2, 4,5)
  K4 <- dplyr::select(star_fusion,Fusion,Gene1,Gene2,base1,chrom1,strand1,base2,chrom2,strand2,spanning.pairs,spanning.reads, GeneId1,GeneId2)
  K4$caller <- "StarFusion"
  TOTAL <- rbind(K1,K2,K3,K4)
  return(TOTAL)
}
TOTAL <- standarise()
write.csv(TOTAL,"1.2practiques/CLEAN_CALLS.csv",row.names = FALSE)
TOTAL <- subselect_from_dups()


selection_bas <- function(one,two){
  if (one == "." && two =="."){
    return(1)
  } 
  if (one == "."){
    return(2)
  }
  if (two == "."){
    return(3)
  }
  return(4)
}
do_lookup <- function(i) {
  if (grepl(",",i)){
    b <- gsub('\\(.*\\).',",",i, perl = TRUE)
    c <- gsub('\\(.*\\)',"",b, perl = TRUE)
    split <- as.data.frame(str_split(c, ","))
    colnames(split)<- "V1"
    total <- c()
    if (grepl("ENSG", split$V1[1])){
      total <- append(total,split$V1[1])
    }
    else{
      tmpS <- xx[xx$symbol == split$V1[1],]
      tmpI <- xxid[xxid$gene_id == tmpS$gene_id,]
      new <- tmpI$ensembl_id
      if (length(new)>1){
        new <- paste(new,collapse=";")
      }
      else if (length(new) == 0){
        new <- "."
      }
      total <- append(total,new)
    }
    if (grepl("ENSG", split$V1[2])){
      total <- append(total,split$V1[2])
    }
    else{
      tmpS <- xx[xx$symbol == split$V1[2],]
      tmpI <- xxid[xxid$gene_id == tmpS$gene_id,]
      new <- tmpI$ensembl_id
      if (length(new)>1){
        new <- paste(new,collapse=";")
      }
      else if (length(new) == 0){
        new <- "."
      }
      total <- append(total,new)
    }
    total <- paste(total, collapse= ",")
    return(total)
  }
  else{
    tmpS <- xx[xx$symbol == as.character(i),]
    tmpI <- xxid[xxid$gene_id == tmpS$gene_id,]
    total <- tmpI$ensembl_id
    if (length(total)>1){
      total <- paste(total,collapse=",")
    }
    else if (length(total) == 0){
      total <- "."
    }
    return(total)
  }
}

chose_id <- function(Genename1,Genename2,GeneId1,GeneId2, range){
  ids <- c()
  if (range == 1){
    Id1 <- do_lookup(Genename1)
    Id2 <- do_lookup(Genename2)
  }
  else if (range == 2){
    Id1 <- do_lookup(Genename1)
    Id2 <- GeneId2
  }
  else{
    Id1 <- GeneId1
    Id2 <- do_lookup(Genename2)
  }
  ids <- append(ids,Id1)
  ids <- append(ids,Id2)
  return(ids)
}

get_ids <- function(Row){
  ids <- c()
  Gene1name <- Row[2]
  Gene2name <- Row[3]
  GeneId1name <- Row[12]
  GeneId2name <- Row[13]
  range <- selection_bas(GeneId1name,GeneId2name)
  if (range == 4){
    ids <- append(ids,GeneId1name)
    ids <- append(ids,GeneId2name)
    return(ids)
  }
  else{
    return(chose_id(Gene1name,Gene2name,GeneId1name,GeneId2name,range))
  }
}


  mapa <-  org.Hs.egSYMBOL
  mapa2 <- org.Hs.egENSEMBL
  mapped_ids <- mappedkeys(mapa2)
  xxid <- as.data.frame(mapa2[mapped_ids])
  mapped_genes <- mappedkeys(mapa)
  # Convert to a list
  xx <- as.data.frame(mapa[mapped_genes])

annotate_ids <- function(){
  for (Row in 1:nrow(TOTAL)){
    ids <- get_ids(TOTAL[Row,])
    TOTAL$GeneId1[Row] <-as.character(ids[1])
    TOTAL$GeneId2[Row] <-as.character(ids[2])
  }
  return(TOTAL)
}



subselect_from_dups <- function(){
 TOTAL <- TOTAL %>% group_by(Fusion,caller)  %>% mutate(Callers = n()) 
 TOTAL <- TOTAL %>% group_by(Fusion) %>% top_n(1,spanning.reads) %>% slice(1)
 return(TOTAL)
}
TOTAL <- annotate_ids()
TOTAL$chrom1[TOTAL$chrom1 == "M"] <- "MT"
TOTAL$chrom2[TOTAL$chrom2 == "M"] <- "MT"
write.csv(TOTAL,"1.2practiques/Standarized_calls.csv" ,row.names = FALSE)
fusion_reads <- read.csv("1.2practiques/BREAK/Results_bam_0.csv")
dataM <- read.csv("1.2practiques/Standarized_calls.csv")
dataM$Proven1 <- fusion_reads$Fused1
dataM$Proven2 <- fusion_reads$Fused2
# sample <- "C11ORF9_RELA.SRR15842441"
# TruePositives = Tp_checking(sample,"edgreen")
# TruePositives
# names(TruePositives) <- "Fusion"
# dataM$Result <- "Novel"
# dataM$Result[dataM$Fusion %in% generalM$Gene] <- "Known"
# dataM$Result[dataM$Fusion %in% TruePositives$Fusion] <- "TP"

pseudogenes <- read.csv("1.2practiques/Pseudos.txt", sep = "\t", header = FALSE)
colnames(pseudogenes)<-  c("Gene1","P1", "Gene2", "P2")
pseudogenes$Fusion <- paste(pseudogenes$Gene1, pseudogenes$Gene2, sep = "::")
Repetitiv <- read.csv("1.2practiques/REPEATED_REG.txt", sep = "\t", header = FALSE)
colnames(Repetitiv)<-  c("Gene1","R1", "Gene2", "R2")
Repetitiv$Fusion <- paste(Repetitiv$Gene1, Repetitiv$Gene2, sep = "::")





