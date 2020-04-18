#Input: Counts RES_DISGENET, RES_KEGG_mart_export, RES_MendelianRanges, RES_OMIM.Entry.Retrieval_monogen3, RES_OMIM.Entry.Retrieval_singlegene3
#These are the output of slippage_finding.py (avaiable on repository)

View(RES_DISGENET)
View(RES_KEGG_mart_export)
View(RES_MendelianRanges)
View(RES_OMIM.Entry.Retrieval_monogen3)
View(RES_OMIM.Entry.Retrieval_singlegene3)

All<- NULL
All$anott<- RES_DISGENET$anott
All$values<- RES_DISGENET$V3 + RES_KEGG_mart_export$V3 +RES_MendelianRanges$V3 + RES_OMIM.Entry.Retrieval_monogen3#V3 + RES_OMIM.Entry.Retrieval_singlegene3$V3


for (i in 1:25){ 
  crom <- All %>% filter(V1== i) 
  CDS= crom %>% filter(anott=="CDS")
  Other=(sum(crom$values) - CDS$values)/sum(crom$values)
  print(i)
  data <- data.frame(
    name=c("CDS","Others") ,  
    value=c((CDS$values)/sum(crom$values),Other)
  )
  if (i==24){
    i="X"
  }
  if (i==25){
    i="Y"
  }
  png(filename = paste(paste("CDS_prop","genomscale",sep = ""),".png",sep = ""),
      width = 572, height = 342, units = "px")
  ##This code represents the core ggplot function for every chromosome graphic, the above part is for the relative proportion of CDS per chromosome, the core ggplot for the ##first plots uses x=annot and with out fill
  print(ggplot(data,aes(x=name,y=value,fill=as.factor(name))) + geom_bar(stat = "identity") 
        +  geom_text(aes(label=value),position=position_dodge(width=0.9), vjust=-0.25)  
        +coord_cartesian(ylim=c(0,1)) 
        +ggtitle("Monogenic proportion - Genomic scale") 
        + labs(x="Annotation",y="Proportion")
        + theme(legend.position = "none"))
  dev.off()
}




some<- p %>% filter(add!=0)
som_all<- some[1:19,]
som_all[20,]<-c("other",20)
##plot for the genomic scale distrbiution of meaning full occurances  
ggplot(som_all,aes(x=class,y=add/sum(add))) + geom_bar(stat = "identity",fill="darkgreen") + theme(axis.text.x = element_text(angle = 0, hjust = 1)) +ggtitle("Genomic elements counts distribution") + labs(x="Annotation",y="Proportion") + coord_flip()


lnc<- All %>% filter(anott=="lnc_RNA")
lnc_val<-sum(lnc$values)
aloth<-sum(All$values)- lnc_val
## with this two elements the lncrna regarding plots


