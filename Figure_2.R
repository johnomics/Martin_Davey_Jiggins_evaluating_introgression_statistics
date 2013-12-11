#!/usr/bin/env Rscript

# Figure_2.R
# Generate the plots (but not the tree diagram) for Figure 2 as a PDF file

# Written for "Evaluating the use of ABBA-BABA statistics to locate introgressed loci"
# by Simon H. Martin, John W. Davey and Chris D. Jiggins
# Simon Martin: shm45@cam.ac.uk
# John Davey:   jd626@cam.ac.uk
# November-December 2013

library(optparse)
suppressMessages(library(ggplot2))
suppressMessages(library(reshape))
library(plyr)
library(stringr)
library(tools)
library(grid)

options<-list(
    make_option(c("-i","--input"),default=NULL,type="character",help="Model CSV file")
)

opt<-parse_args(OptionParser(option_list=options))
if (is.null(opt$input)) stop("Please supply an input file with -i")
if (!file.exists(opt$input)) stop(opt$input, "does not exist! Please specify a valid model CSV file with -i")

partition.x<-function(model,stat,threshold=0.1,win) {
    outnum<-threshold*win
    Partition<-sapply(rank(-model[,stat]), function(x) ifelse(x<=outnum,"_Outlier","_Background"))
    Partition<-paste0(stat,Partition)
    model<-cbind(model,Partition)
    names(model)[names(model)=="Partition"]<-paste0(stat,"_Partition")
    model
}

model<-read.csv(opt$input)

pars<-data.frame(str_match(opt$input,"Alternate_t123-(.+)_t23-(.+)\\.Background_t123-(.+)_t21-(.+)\\.(.+)\\.csv"),stringsAsFactors=FALSE)
names(pars)<-c("File","Alternate_t123","Alternate_t23","Background_t123","Background_t21","Windows")
pars$Windows<-as.numeric(pars$Windows)

model$Model<-factor(str_match(model$Model,"(.+)_t123")[,2],levels=c("Background","Alternate"))

model<-partition.x(model,"D",0.1,pars$Windows)
model<-partition.x(model,"mfD0",0.1,pars$Windows)

model$ModelD0<-ifelse(model$D>=0,as.character(model$Model),"Unused")

D.threshold<-min(model[model$D_Partition=="D_Outlier",]$D)
f.threshold<-min(model[model$mfD0_Partition=="mfD0_Outlier",]$mfD0)

pdf("Figure_2_plots.pdf",8,7)
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,2)))

treecols<-c("Background"="#ad5c08","Alternate"="#00ba38","Unused"="white")

d.plot<-ggplot(model,aes(D,colour=Model,fill=ModelD0))+
        geom_histogram(binwidth=0.07475,position="dodge")+
        geom_vline(aes(xintercept=D.threshold),linetype="longdash")+
        xlim(c(-1.03,1.03))+
        ylab("Count")+
        scale_fill_manual(values=treecols)+
        scale_colour_manual(values=treecols)+
        guides(fill=FALSE,colour=FALSE)+
        theme_bw(base_family="Helvetica")+
        theme(axis.title=element_text(face="bold"))

f.plot<-ggplot(model,aes(mfD0,fill=Model))+
        geom_histogram(binwidth=0.072,position="dodge")+
        geom_vline(aes(xintercept=f.threshold),linetype="longdash")+
        xlim(c(0,1))+
        xlab("f")+
        ylab("Count")+
        scale_fill_manual(values=treecols)+
        guides(fill=FALSE)+
        theme_bw(base_family="Helvetica")+
        theme(axis.title.x=element_text(face="bold.italic"),axis.title.y=element_text(face="bold"))

print(d.plot,vp=viewport(layout.pos.row=2,layout.pos.col=c(1,2)))
print(f.plot,vp=viewport(layout.pos.row=1,layout.pos.col=2))

dev.off()
q()


# Preserve order by Model, then sort by D
model.block<-model[c(1:100,1001:1900),]
model.block<-model.block[order(model.block$D,na.last=FALSE),]
model.block<-melt(model.block[,c("D","mfD0","Model")],id=c("Model"))
model.block$Windows<-1:1000
model.block$variable<-revalue(model.block$variable,c(mfD0="f"))


model.plot<-model[order(model$Model,model$D,na.last=FALSE),]
model.plot<-melt(model.plot[,c("D","mfD0","Model")],id=c("Model"))
model.plot$Windows<-1:10000
model.plot$variable<-revalue(model.plot$variable,c(mfD0="f"))

ggplot(model.plot,aes(value,fill=Model))+
      geom_histogram(binwidth=0.07,position="dodge")+
      theme_bw()+
      geom_vline(aes(xintercept=thresholds),data=data.frame(variable=c("D","f"),thresholds=c(D.threshold,f.threshold)))+
      facet_grid(variable~.,scales="free_y")

ggplot(model.plot,aes(Windows,colour=Model),environment=environment())+
        geom_segment(aes(x=Windows,xend=Windows,y=0,yend=value),environment=environment())+
        theme_bw()+
        geom_hline(aes(yintercept=thresholds),data=data.frame(variable=c("D","f"),thresholds=c(D.threshold,f.threshold)))+
        facet_grid(variable~.,scales="free_y")+
        theme(axis.title.y=element_blank())+aes(ymax=1)

ggplot(model.block,aes(Windows,colour=Model),environment=environment())+
        geom_segment(aes(x=Windows,xend=Windows,y=0,yend=value),environment=environment())+
        theme_bw()+
        geom_hline(aes(yintercept=thresholds),data=data.frame(variable=c("D","f"),thresholds=c(D.threshold,f.threshold)))+
        facet_grid(variable~.,scales="free_y")+
        theme(axis.title.y=element_blank())+aes(ymax=1)


dev.off()
