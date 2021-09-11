library(pheatmap)
library(ggplot2)
library(S4Vectors)

#TCR use baise
setwd("E:/OSCC/integrate/3_5/FastMNN/T/T.del.double.dead/CD4T.without.RPL")
load("CD4T.RData")
CD4T.meta<-CD4T@meta.data
write.csv(CD4T.meta,file = "CD4T.meta.csv")
#TRA
TRAV<-CD4T.meta$alpha.V
TRAV<-unlist(strsplit(TRAV,","))
TRAV<-TRAV[!duplicated(TRAV)]
TRAV<-TRAV[!is.na(TRAV)]
TRAV<-TRAV[-which(TRAV=="NA")]
TRAV.df<-as.data.frame(array(0,c(8,length(TRAV))))
colnames(TRAV.df)<-TRAV
row.names(TRAV.df)<-names(table(CD4T@meta.data$Celltype))

TRAJ<-CD4T.meta$alpha.J
TRAJ<-unlist(strsplit(TRAJ,","))
TRAJ<-TRAJ[!duplicated(TRAJ)]
TRAJ<-TRAJ[!is.na(TRAJ)]
TRAJ<-TRAJ[-which(TRAJ=="NA")]
TRAJ.df<-as.data.frame(array(0,c(8,length(TRAJ))))
colnames(TRAJ.df)<-TRAJ
row.names(TRAJ.df)<-names(table(CD4T@meta.data$Celltype))

TRAC<-CD4T.meta$alpha.C
TRAC<-unlist(strsplit(TRAC,","))
TRAC<-TRAC[!duplicated(TRAC)]
TRAC<-TRAC[!is.na(TRAC)]
TRAC<-TRAC[-which(TRAC=="NA")]
TRAC.df<-as.data.frame(array(0,c(8,length(TRAC))))
colnames(TRAC.df)<-TRAC
row.names(TRAC.df)<-names(table(CD4T@meta.data$Celltype))

#TRB
TRBV<-CD4T.meta$beta.V
TRBV<-unlist(strsplit(TRBV,","))
TRBV<-TRBV[!duplicated(TRBV)]
TRBV<-TRBV[!is.na(TRBV)]
TRBV<-TRBV[-which(TRBV=="NA")]
TRBV.df<-as.data.frame(array(0,c(8,length(TRBV))))
colnames(TRBV.df)<-TRBV
row.names(TRBV.df)<-names(table(CD4T@meta.data$Celltype))

TRBD<-CD4T.meta$beta.D
TRBD<-unlist(strsplit(TRBD,","))
TRBD<-TRBD[!duplicated(TRBD)]
TRBD<-TRBD[!is.na(TRBD)]
TRBD<-TRBD[-which(TRBD=="NA")]
TRBD.df<-as.data.frame(array(0,c(8,length(TRBD))))
colnames(TRBD.df)<-TRBD
row.names(TRBD.df)<-names(table(CD4T@meta.data$Celltype))

TRBJ<-CD4T.meta$beta.J
TRBJ<-unlist(strsplit(TRBJ,","))
TRBJ<-TRBJ[!duplicated(TRBJ)]
TRBJ<-TRBJ[!is.na(TRBJ)]
TRBJ<-TRBJ[-which(TRBJ=="NA")]
TRBJ.df<-as.data.frame(array(0,c(8,length(TRBJ))))
colnames(TRBJ.df)<-TRBJ
row.names(TRBJ.df)<-names(table(CD4T@meta.data$Celltype))

TRBC<-CD4T.meta$beta.C
TRBC<-unlist(strsplit(TRBC,","))
TRBC<-TRBC[!duplicated(TRBC)]
TRBC<-TRBC[!is.na(TRBC)]
TRBC<-TRBC[-which(TRBC=="NA")]
TRBC.df<-as.data.frame(array(0,c(8,length(TRBC))))
colnames(TRBC.df)<-TRBC
row.names(TRBC.df)<-names(table(CD4T@meta.data$Celltype))

for(cluster in names(table(CD4T@meta.data$Celltype))){
  #cluster="CD4T_tra_MKI67"
  TRAV=CD4T.meta[which(CD4T.meta$Celltype==cluster),"alpha.V"]
  TRAV<-unlist(strsplit(TRAV,","))
  TRAV<-TRAV[!is.na(TRAV)]
  if(!isEmpty(which(TRAV=="NA"))){
    TRAV<-TRAV[-which(TRAV=="NA")]
  }
  count=table(TRAV)
  for (col in colnames(TRAV.df)[1:ncol(TRAV.df)]) {
    TRAV.df[as.character(cluster),col]=count[col]
  }
  TRAV.df[is.na(TRAV.df)]<-0
  TRAV.df[as.character(cluster),]=TRAV.df[as.character(cluster),]/sum(TRAV.df[as.character(cluster),])
  
  TRAJ=CD4T.meta[which(CD4T.meta$Celltype==cluster),"alpha.J"]
  TRAJ<-unlist(strsplit(TRAJ,","))
  TRAJ<-TRAJ[!is.na(TRAJ)]
  if(!isEmpty(which(TRAJ=="NA"))){
    TRAJ<-TRAJ[-which(TRAJ=="NA")]
  }
  count=table(TRAJ)
  for (col in colnames(TRAJ.df)[1:ncol(TRAJ.df)]) {
    TRAJ.df[as.character(cluster),col]=count[col]
  }
  TRAJ.df[is.na(TRAJ.df)]<-0
  TRAJ.df[as.character(cluster),]=TRAJ.df[as.character(cluster),]/sum(TRAJ.df[as.character(cluster),])
  
  TRAC=CD4T.meta[which(CD4T.meta$Celltype==cluster),"alpha.C"]
  TRAC<-unlist(strsplit(TRAC,","))
  TRAC<-TRAC[!is.na(TRAC)]
  if(!isEmpty(which(TRAC=="NA"))){
    TRAC<-TRAC[-which(TRAC=="NA")]
  }
  count=table(TRAC)
  for (col in colnames(TRAC.df)[1:ncol(TRAC.df)]) {
    TRAC.df[as.character(cluster),col]=count[col]
  }
  TRAC.df[is.na(TRAC.df)]<-0
  TRAC.df[as.character(cluster),]=TRAC.df[as.character(cluster),]/sum(TRAC.df[as.character(cluster),])
  
  TRBV=CD4T.meta[which(CD4T.meta$Celltype==cluster),"beta.V"]
  TRBV<-unlist(strsplit(TRBV,","))
  TRBV<-TRBV[!is.na(TRBV)]
  if(!isEmpty(which(TRBV=="NA"))){
    TRBV<-TRBV[-which(TRBV=="NA")]
  }
  count=table(TRBV)
  for (col in colnames(TRBV.df)[1:ncol(TRBV.df)]) {
    TRBV.df[as.character(cluster),col]=count[col]
  }
  TRBV.df[is.na(TRBV.df)]<-0
  TRBV.df[as.character(cluster),]=TRBV.df[as.character(cluster),]/sum(TRBV.df[as.character(cluster),])
  
  TRBD=CD4T.meta[which(CD4T.meta$Celltype==cluster),"beta.D"]
  TRBD<-unlist(strsplit(TRBD,","))
  TRBD<-TRBD[!is.na(TRBD)]
  if(!isEmpty(which(TRBD=="NA"))){
    TRBD<-TRBD[-which(TRBD=="NA")]
  }
  count=table(TRBD)
  for (col in colnames(TRBD.df)[1:ncol(TRBD.df)]) {
    TRBD.df[as.character(cluster),col]=count[col]
  }
  TRBD.df[is.na(TRBD.df)]<-0
  TRBD.df[as.character(cluster),]=TRBD.df[as.character(cluster),]/sum(TRBD.df[as.character(cluster),])
  
  TRBJ=CD4T.meta[which(CD4T.meta$Celltype==cluster),"beta.J"]
  TRBJ<-unlist(strsplit(TRBJ,","))
  TRBJ<-TRBJ[!is.na(TRBJ)]
  if(!isEmpty(which(TRBJ=="NA"))){
    TRBJ<-TRBJ[-which(TRBJ=="NA")]
  }
  count=table(TRBJ)
  for (col in colnames(TRBJ.df)[1:ncol(TRBJ.df)]) {
    TRBJ.df[as.character(cluster),col]=count[col]
  }
  TRBJ.df[is.na(TRBJ.df)]<-0
  TRBJ.df[as.character(cluster),]=TRBJ.df[as.character(cluster),]/sum(TRBJ.df[as.character(cluster),])
  
  TRBC=CD4T.meta[which(CD4T.meta$Celltype==cluster),"beta.C"]
  TRBC<-unlist(strsplit(TRBC,","))
  TRBC<-TRBC[!is.na(TRBC)]
  if(!isEmpty(which(TRBC=="NA"))){
    TRBC<-TRBC[-which(TRBC=="NA")]
  }
  count=table(TRBC)
  for (col in colnames(TRBC.df)[1:ncol(TRBC.df)]) {
    TRBC.df[as.character(cluster),col]=count[col]
  }
  TRBC.df[is.na(TRBC.df)]<-0
  TRBC.df[as.character(cluster),]=TRBC.df[as.character(cluster),]/sum(TRBC.df[as.character(cluster),])
}

TRAV.df<-TRAV.df[,c("TRAV1-1","TRAV1-2","TRAV2","TRAV3","TRAV4","TRAV5","TRAV6",
                    "TRAV8-1", "TRAV8-2","TRAV8-3","TRAV8-4" ,"TRAV8-6",
                    "TRAV9-2","TRAV10" ,"TRAV12-1","TRAV12-2","TRAV12-3","TRAV13-1",
                    "TRAV13-2","TRAV14/DV4","TRAV16","TRAV17","TRAV18","TRAV19" ,
                    "TRAV20","TRAV21" , "TRAV22" ,"TRAV23/DV6","TRAV24","TRAV25","TRAV26-1",
                    "TRAV26-2","TRAV27" ,"TRAV29/DV5","TRAV30","TRAV34","TRAV35" ,
                    "TRAV36/DV7","TRAV38-1", "TRAV38-2/DV8","TRAV39","TRAV40","TRAV41" )]

TRAJ.df<-TRAJ.df[,c("TRAJ3", "TRAJ4", "TRAJ5", "TRAJ6", "TRAJ7",  "TRAJ8",  "TRAJ9" ,  "TRAJ10", "TRAJ11", "TRAJ12", "TRAJ13", 
                    "TRAJ15", "TRAJ16", "TRAJ17", "TRAJ18", "TRAJ20", "TRAJ8", "TRAJ22", "TRAJ23",
                    "TRAJ24", "TRAJ26", "TRAJ27", "TRAJ28", "TRAJ29", "TRAJ30", "TRAJ31", "TRAJ32", "TRAJ33", "TRAJ34", "TRAJ35",
                    "TRAJ36", "TRAJ37", "TRAJ38", "TRAJ39", "TRAJ40", "TRAJ41", "TRAJ42", "TRAJ43", "TRAJ44", "TRAJ45", "TRAJ46",
                    "TRAJ47", "TRAJ48", "TRAJ49", "TRAJ50", "TRAJ52", "TRAJ53", "TRAJ54", "TRAJ56", "TRAJ57", "TRAJ58", "TRAJ61")]

TRBV.df<-TRBV.df[,c("TRBV2","TRBV3-1","TRBV4-1","TRBV4-2","TRBV4-3","TRBV5-1","TRBV5-3","TRBV5-4","TRBV5-5",
                    "TRBV5-6","TRBV5-8","TRBV6-1","TRBV6-2","TRBV6-4","TRBV6-5","TRBV6-6", 
                    "TRBV6-7","TRBV7-2","TRBV7-3","TRBV7-4","TRBV7-6","TRBV7-7",
                    "TRBV7-8","TRBV7-9","TRBV9","TRBV10-1","TRBV10-2","TRBV10-3","TRBV11-1",
                    "TRBV11-2","TRBV11-3","TRBV12-3","TRBV12-4","TRBV12-5","TRBV13","TRBV14",
                    "TRBV15","TRBV16","TRBV18","TRBV19","TRBV20-1","TRBV21-1","TRBV23-1","TRBV24-1",
                    "TRBV25-1","TRBV27","TRBV28","TRBV29-1","TRBV30")]

TRBJ.df<-TRBJ.df[,c("TRBJ1-1",  "TRBJ1-2",  "TRBJ1-3",  "TRBJ1-4",  "TRBJ1-5",  "TRBJ1-6",  "TRBJ2-1",  
                    "TRBJ2-2",  "TRBJ2-2P", "TRBJ2-3",  "TRBJ2-4",  "TRBJ2-5",  "TRBJ2-6",  "TRBJ2-7" )]
TRBC.df<-TRBC.df[,c("TRBC1","TRBC2")]

#CD4

pdf(file="pheatmap.TRAV.pdf", width = 14, height = 6)  
pheatmap(TRAV.df,treeheight_row=0,treeheight_col=0,angle_col=45,cluster_rows =F,
         color = colorRampPalette(c("#4F8ABC", "white","#D22D37"))(50),cluster_cols =F,
         border_color = "black",cellwidth = 15, cellheight =15)
dev.off()

pdf(file="bar.TRAV23.DV6_TRAJ23.pdf", width = 5, height = 3)  
ggplot(CD4T.meta[which(CD4T.meta$alpha.V=="TRAV23/DV6" & CD4T.meta$Celltype=="CD4T_exh_LAG3"),],aes(x=clonotype.10x,fill=sample))+
  geom_bar(stat="count")+
  geom_text(stat='count', aes(label=..count..), hjust=0)+
  scale_fill_manual(values=c("#BC8F8F","#20B2AA","#FF6A6A"))+
  theme(axis.text.x = element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  coord_flip()+
  labs(x="",y="alpha.V=TRAV23/DV6")
ggplot(CD4T.meta[which(CD4T.meta$alpha.J=="TRAJ23" & CD4T.meta$Celltype=="CD4T_exh_LAG3"),],aes(x=clonotype.10x,fill=sample))+
  geom_bar(stat="count")+
  geom_text(stat='count', aes(label=..count..), hjust=0)+
  scale_fill_manual(values=c("#20B2AA"))+
  theme(axis.text.x = element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  coord_flip()+
  labs(x="",y="alpha.J=TRAJ23")
dev.off()

pdf(file="pheatmap.TRAJ.pdf", width = 16, height = 6)  
pheatmap(TRAJ.df,treeheight_row=0,treeheight_col=0,angle_col=45,cluster_rows =F,
         color = colorRampPalette(c("#4F8ABC", "white","#D22D37"))(50),cluster_cols =F,
         border_color = "black",cellwidth = 15, cellheight =15)
dev.off()

pdf(file="pheatmap.TRBV.pdf", width = 16, height = 6)  
pheatmap(TRBV.df,treeheight_row=0,treeheight_col=0,angle_col=45,cluster_rows =F,
         color = colorRampPalette(c("#4F8ABC", "white","#D22D37"))(50),cluster_cols =F,
         border_color = "black",cellwidth = 15, cellheight =15)
dev.off()

pdf(file="pheatmap.TRBD.pdf", width = 8, height = 6)  
pheatmap(TRBD.df,treeheight_row=0,treeheight_col=0,angle_col=45,cluster_rows =F,
         color = colorRampPalette(c("#4F8ABC", "white","#D22D37"))(50),cluster_cols =F,
         border_color = "black",cellwidth = 15, cellheight =15)
dev.off()

pdf(file="pheatmap.TRBJ.pdf", width = 8, height = 6)  
pheatmap(TRBJ.df,treeheight_row=0,treeheight_col=0,angle_col=45,cluster_rows =F,
         color = colorRampPalette(c("#4F8ABC", "white","#D22D37"))(50),cluster_cols =F,
         border_color = "black",cellwidth = 15, cellheight =15)
dev.off()

pdf(file="pheatmap.TRBC.pdf", width = 8, height = 6)  
pheatmap(TRBC.df,treeheight_row=0,treeheight_col=0,angle_col=45,cluster_rows =F,
         color = colorRampPalette(c("#4F8ABC", "white","#D22D37"))(50),cluster_cols =F,
         border_color = "black",cellwidth = 15, cellheight =15)
dev.off()

pdf(file="bar.TRBV6-1.TRBJ1-2.pdf", width = 5, height = 3)  
ggplot(CD4T.meta[which(CD4T.meta$beta.V=="TRBV6-1" & CD4T.meta$Celltype=="CD4T_exh_LAG3"),],aes(x=clonotype.10x,fill=sample))+
  geom_bar(stat="count")+
  geom_text(stat='count', aes(label=..count..), hjust=0)+
  scale_fill_manual(values=c("#BC8F8F","#20B2AA","#FF6A6A"))+
  theme(axis.text.x = element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  coord_flip()+
  labs(x="",y="beta.V=TRBV6-1")
ggplot(CD4T.meta[which(CD4T.meta$beta.J=="TRBJ1-2" & CD4T.meta$Celltype=="CD4T_exh_LAG3"),],aes(x=clonotype.10x,fill=sample))+
  geom_bar(stat="count")+
  geom_text(stat='count', aes(label=..count..), hjust=0)+
  scale_fill_manual(values=c("#BC8F8F","#20B2AA","#FF6A6A"))+
  theme(axis.text.x = element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  coord_flip()+
  labs(x="",y="beta.J=TRBJ1-2")
dev.off()
