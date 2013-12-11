#!/usr/bin/env Rscript

# Figure_S1_S2.R
# Generate the plots for Figures S1 and S2 as PDF files

# Written for "Evaluating statistics for the identification of introgressed loci"
# by Simon H. Martin, John W. Davey and Chris D. Jiggins
# Simon Martin: shm45@cam.ac.uk
# John Davey:   jd626@cam.ac.uk
# November-December 2013


library(optparse)
suppressMessages(library(ggplot2))
library(plyr)

options<-list(
    make_option(c("-i","--input"),default=NULL,type="character",help="Input file")
)

opt<-parse_args(OptionParser(option_list=options))
if (is.null(opt$input)) stop("Please supply an input file with -i")

make.grid.plot<-function(models,stat) {
    modelset<-subset(models,Stat==stat)
    modelset$Alternate_t123<-factor(modelset$Alternate_t123,levels=seq(2,0.2,-0.2))

    modelset$ModelType<-revalue(modelset$ModelType,c("Ancestral structure"="Ancestral Structure","gene flow P2-P3"="Gene Flow P2-P3","gene flow P3-P2"="Gene Flow P3-P2"))
    modelset$Stat<-revalue(modelset$Stat,c(mfD0="f"))

    names(modelset)[1]<-"Bgt123"
    names(modelset)[3]<-"Alt123"
    ggplot(modelset,
    	   aes(Background_t21,Alternate_t23,size=OutlierAlternatePC,colour=ModelType)
    	  )+
    	   geom_point()+
    	   facet_grid(Alt123~Bgt123,labeller=label_both)+
    	   scale_size(paste("% predicted by",unique(modelset$Stat)),limits=c(0,100),breaks=seq(0,100,20),range=c(1,9))+
    	   scale_x_continuous(expression(paste(Bgt[21])),breaks=seq(0,2,0.2))+
    	   scale_y_continuous(expression(paste(Alt[23])),breaks=seq(0,2,0.2))+
    	   theme_bw(base_size=10)
}


models<-read.delim(opt$input)

pdf("Figure_S1.pdf",21,21)
make.grid.plot(models,"D")
dev.off()

pdf("Figure_S2.pdf",21,21)
make.grid.plot(models,"mfD0")
dev.off()

warninglist<-warnings()
if (!(is.null(warninglist))) warninglist