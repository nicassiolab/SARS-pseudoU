library(tidyverse)
library(dplyr)
library(magic)

ROOTDIR="/Volumes/scratch/FN/TL/cugolini/cov/analysis/miseq/map_to_genome"

# Function the returns full path from basedir
bdp <- function(relpath){
  return(paste0(ROOTDIR,"/",relpath))
}

bedfile <-
  read.table(
    file = "Aligned.out.bed",
    col.names = c(
      "chr",
      "start",
      "end",
      "name",
      "score",
      "strand",
      "thickStart",
      "thickEnd",
      "itemRgb",
      "blockCount",
      "blockSizes",
      "blockStarts"
    )
  )

two_ex <- bedfile %>%
  subset(blockCount==2)%>%
  separate(blockStarts,into = c("start1","start2"),sep = ",")%>%
  separate(blockSizes,into = c("size1","size2"),sep = ",")%>%
  mutate(start1=(as.numeric(start1)+as.numeric(start)))%>%
  mutate(start2=(as.numeric(start2)+as.numeric(start)))%>%
  mutate(firstblockend=as.numeric(start)+as.numeric(size1))

breaks <- seq(0, 30000, by = 1000)
intervals <- data.frame(breaks,shift(breaks,i=-1))
intervals <- head(intervals,-1)
colnames(intervals)<-c("start","end")

df <- data.frame(
                 juncstart=character(),
                 freq=double(),
                 juncstart=character(),
                 stringsAsFactors=FALSE)

for(i in 1:nrow(intervals)) {
  appo <- subset(two_ex,firstblockend>=intervals$start[i] & firstblockend<intervals$end[i])
  juncend <-as.data.frame(table(cut(appo$start2, breaks = breaks)))
  pippo <- cbind(juncend,juncend$Var1[i])
  colnames(pippo)<- c("juncend","freq","juncstart")
  df<-rbind(df,pippo)
  print(length(appo$start2))
}


p <- df %>%
  {
    ggplot(., aes(x=juncend, y=juncstart, fill=log10(freq+1))) +
      geom_tile() +
      scale_fill_distiller(name="log10(Frequency)", palette="GnBu", direction=1) +
      theme_bw(5) + vlab 
  }
p



breaks <- seq(30, 100, by = 10)
breaks2 <- seq(28980, 29800, by = 20)
breaks2_middlepoint<-seq(28990,29790,by=20)
intervals <- data.frame(breaks,shift(breaks,i=-1))
intervals <- head(intervals,-1)
colnames(intervals)<-c("start","end")

df <- data.frame(
  juncstart=character(),
  freq=double(),
  juncstart=character(),
  stringsAsFactors=FALSE)

for(i in 1:nrow(intervals)) {
  appo <- subset(two_ex,firstblockend>=intervals$start[i] & firstblockend<intervals$end[i])
  juncend <-as.data.frame(table(cut(appo$start2, breaks = breaks2)))
  pippo <- cbind(juncend,paste0("[",intervals$start[i],",",intervals$end[i],")"),as.character(breaks2_middlepoint))
  colnames(pippo)<- c("juncend","freq","juncstart","juncend_middlepoint")
  df<-rbind(df,pippo)
  print(length(appo$start2))
}


p <- df %>%
  {
    ggplot(., aes(x=juncend_middlepoint, y=juncstart, fill=log10(freq+1))) +
      geom_tile() +
      scale_fill_distiller(name="Log10(Frequency)", palette="GnBu", direction=1) +
      theme_bw(30) + vlab 
  }
p
ggsave(p,file=bdp("junctions_heatmap.eps"), width = 35,height = 20)


orf9d_extr<- subset(two_ex,firstblockend>=30 & firstblockend<100) %>%
  subset(start2>=29100 & start2<29180)

orf9d_extr<- subset(bedfile,bedfile$name %in% orf9d_extr$name)
write.table(orf9d_extr, sep = "\t",quote=F,row.names = F,col.names = F, file=bdp("orf9d.bed"))


orf10_int_extr<- subset(two_ex,firstblockend>=30 & firstblockend<100) %>%
  subset(start2>=29610 & start2<29630)
orf10_int_extr<-head(orf10_int_extr,150) #extract first 150 reads as an example for the track
orf10_est<- subset(bedfile,bedfile$name %in% orf10_int_extr$name)
write.table(orf10_int, sep = "\t",quote=F,row.names = F,col.names = F, file=bdp("orf10_int.bed"))


orf10_est_extr<- subset(two_ex,firstblockend>=30 & firstblockend<100) %>%
  subset(start2>=29500 & start2<29540)
orf10_est_extr<-head(orf10_est_extr,150)
orf10_est<- subset(bedfile,bedfile$name %in% orf10_est_extr$name)
write.table(orf10_est, sep = "\t",quote=F,row.names = F,col.names = F, file=bdp("orf10_est.bed"))

