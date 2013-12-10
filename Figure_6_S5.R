#!/usr/bin/env Rscript

suppressMessages(library(ggplot2))
suppressMessages(library(reshape))
library(plyr)
library(grid)

options(width=1000)

m.all.equal<-Vectorize(function(x, y) {isTRUE(all.equal(x, y))})

get.sig.val<-function(pval) {
    val<-""
    if (pval >= 0.05) val<-"#"
    if (pval <0.05) val<-"*"
    if (pval<0.01) val<-"**"
    if (pval<0.001) val<-"***"
    val        
}

get.sig.text<-function(df,mt) {
    sig.text<-ddply(df, .(Stat,Comparison), function(x) {
        val<-""
        if (mt == "gene flow P2-P3") {
            if (x$Comparison[1] == "P2P3") val<-get.sig.val(unique(x$p.background.higher))
            if (x$Comparison[1] == "P1P3") val<-get.sig.val(unique(x$p.background.higher))
            if (x$Comparison[1] == "P1P2") val<-get.sig.val(unique(x$p.two.tailed))
        }
        if (mt == "Ancestral structure") {
            if (x$Comparison[1] == "P2P3") val<-get.sig.val(unique(x$p.two.tailed))
            if (x$Comparison[1] == "P1P3") val<-get.sig.val(unique(x$p.background.lower))
            if (x$Comparison[1] == "P1P2") val<-get.sig.val(unique(x$p.background.lower))
        }
        val
    })
    names(sig.text)<-c("Stat","Comparison","SigText")
    sig.text
}

skplot<-function(df,mt) {
    
    max.val<-max(df$Mean+df$SE)
    # Add significance line positions
    df<-cbind(df,
              SigXpos1=rep(1                             ,nrow(df)  ),
              SigYpos1=rep(c(max.val+0.005,max.val+0.01) ,nrow(df)/2),
              SigXpos2=rep(c(1,2)                        ,nrow(df)/2),
              SigYpos2=rep(max.val+0.01                  ,nrow(df)  ),
              SigXpos3=rep(2                             ,nrow(df)  ),
              SigYpos3=rep(c(max.val+0.01,max.val+0.005) ,nrow(df)/2),
              max.val=max.val
        )

    # Get significance symbols (*, # etc)
    df<-merge(df,get.sig.text(df,mt),by=c("Stat","Comparison"))

    df$Comparison<-revalue(df$Comparison, c(P1P2="P1-P2",P1P3="P1-P3",P2P3="P2-P3"))

    ggplot(df,aes(Partition,Mean))+
                  facet_grid(Stat ~ Comparison)+
                  geom_bar(stat="identity",colour="black",fill=rep(c("white","grey")))+
                  geom_errorbar(aes(ymax=Mean+SE,ymin=Mean-SE,width=0.25))+
                  geom_path(aes(SigXpos1,SigYpos1))+
                  geom_path(aes(SigXpos2,SigYpos2))+
                  geom_path(aes(SigXpos3,SigYpos3))+
                  geom_text(aes(label=SigText,x=1.5,y=max.val+0.015),size=4)+
                  ylim(c(0,max.val+0.025))+
                  ylab(expression(Mean~italic(d[XY])))+
                  theme_bw(base_size=10,base_family="Helvetica")+
                  theme(axis.text.x=element_text(size=6),
                        strip.text=element_text(size=8),
                        plot.margin=unit(c(0,0,0,0),"mm"))
}

get.prediction<-function(mt) {

    Mean<-NULL

    if (mt=="Ancestral structure") {
        Mean<-c(0.045,0.045,0.05,0.045,0.05,0.035)
        p.two.tailed<-c(rep(0.1,2),rep(0.0001,2),rep(0.0001,2))
        p.background.higher<-c(rep(0.1,2),rep(0.1,2),rep(0.1,2))
        p.background.lower<-c(rep(0.1,2),rep(0.0001,2),rep(0.0001,2))
    } else if (mt=="gene flow P2-P3") {
        Mean<-c(0.025,0.05,0.035,0.05,0.035,0.035)
        p.two.tailed<-c(rep(0.0001,2),rep(0.0001,2),rep(0.1,2))
        p.background.higher<-c(rep(0.0001,2),rep(0.0001,2),rep(0.1,2))
        p.background.lower<-c(rep(0.1,2),rep(0.1,2),rep(0.1,2))
    }
    if (is.null(Mean)) stop("Can't define Prediction")

    data.frame(Stat=rep("Prediction",6),
               Partition=rep(c("Outlier","Background"),3),
               Comparison=factor(c("P2P3","P2P3","P1P3","P1P3","P1P2","P1P2"),levels=c("P2P3","P1P3","P1P2")),
               Mean=Mean,
               SE=rep(0,6),
               p.two.tailed=p.two.tailed,
               p.background.higher=p.background.higher,
               p.background.lower=p.background.lower)
}

get.dxy.df<-function(dxy.stats) {
    dxy.df<-dxy.stats[grep("dxy",dxy.stats$Variable),c("Model","Variable","Mean","SE","w.p.two.tailed","w.p.background.higher","w.p.background.lower")]
    dxy.df<-cbind(dxy.df,colsplit(dxy.df$Variable,split="_",names=c("Comparison","dxy")),colsplit(dxy.df$Model,split="_",names=c("Stat","Partition")))
    dxy.df$Comparison<-factor(dxy.df$Comparison,levels=c("P2P3","P1P3","P1P2"))
    dxy.df$Partition<-factor(dxy.df$Partition,levels=c("Outlier","Background"))
    dxy.df<-rename(dxy.df,c("w.p.two.tailed"="p.two.tailed","w.p.background.higher"="p.background.higher","w.p.background.lower"="p.background.lower"))
    dxy.df[,c("Stat","Partition","Comparison","Mean","SE","p.two.tailed","p.background.higher","p.background.lower")]
}

get.figure5.data<-function(models,bt123,bt21,at123,at23,vpcol) {
    ms<-subset(models,m.all.equal(Background_t123,bt123) & m.all.equal(Background_t21,bt21) & m.all.equal(Alternate_t123,at123) & m.all.equal(Alternate_t23,at23))
    if (nrow(ms) == 0) stop("No model for these parameters")
    mt<-unique(ms$ModelType)
    Prediction<-get.prediction(mt)

    ms.dxy<-get.dxy.df(ms)
    ms.dxy<-rbind(ms.dxy, Prediction)
    ms.dxy$Stat<-factor(ms.dxy$Stat,levels=c("Prediction","Real","D","mfD0"))
    ms.dxy$Stat<-revalue(ms.dxy$Stat, c(Real="Simulation",mfD0="f"))
    ms.dxy$Partition<-revalue(ms.dxy$Partition,c(Outlier="Alternate"))
    print(skplot(ms.dxy[ms.dxy$Stat %in% c("Prediction","Simulation"),],mt),vp=viewport(layout.pos.row=1,layout.pos.col=vpcol))
    ms.dxy$Partition<-revalue(ms.dxy$Partition,c(Alternate="Outlier",Background="Non-Outlier"))
    print(skplot(ms.dxy[ms.dxy$Stat %in% c("D","f"),],mt),vp=viewport(layout.pos.row=2,layout.pos.col=vpcol))
}


# Figure 6 - simulated data
models<-read.delim("model_files_win10000_size5000.dxy.summary.tsv",stringsAsFactors=FALSE)

pdf("Figure_6.pdf",width=6.75,height=6.75)
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,2)))
#get.figure5.data(models,bt123=0.4,bt21=0.2,at123=1.8,at23=0.4,vpcol=1)
#get.figure5.data(models,bt123=2.0,bt21=1.2,at123=1.2,at23=0.4,vpcol=1) # Gene flow P2-P3
get.figure5.data(models,bt123=2.0,bt21=1.2,at123=1.2,at23=0.2,vpcol=1)
get.figure5.data(models,bt123=1.8,bt21=0.8,at123=2,at23=1.8,vpcol=2) # Ancestral structure

grid.text("A",x=0.025,y=0.99,just=c("left","top"),vp=viewport(layout.pos.row=1,layout.pos.col=1),gp=gpar(fontface="bold",fontsize=14))
grid.text("B",x=0.025,y=0.99,just=c("left","top"),vp=viewport(layout.pos.row=2,layout.pos.col=1),gp=gpar(fontface="bold",fontsize=14))
grid.text("C",x=0.025,y=0.99,just=c("left","top"),vp=viewport(layout.pos.row=1,layout.pos.col=2),gp=gpar(fontface="bold",fontsize=14))
grid.text("D",x=0.025,y=0.99,just=c("left","top"),vp=viewport(layout.pos.row=2,layout.pos.col=2),gp=gpar(fontface="bold",fontsize=14))
dev.off()

# Figure S5 - real data
real<-read.delim("Heliconius_genome_windows.dxy.summary.tsv",stringsAsFactors=FALSE)

real.dxy<-get.dxy.df(real)
real.dxy$Stat<-factor(real.dxy$Stat,levels=c("D","mfD0"))
real.dxy$Stat<-revalue(real.dxy$Stat, c(mfD0="f"))
real.dxy$Partition<-revalue(real.dxy$Partition,c(Background="Non-Outlier"))
pdf("Figure_S5.pdf",3.25,3.25)
skplot(real.dxy,"gene flow P2-P3")
dev.off()

warninglist<-warnings()
if (!(is.null(warninglist))) warninglist