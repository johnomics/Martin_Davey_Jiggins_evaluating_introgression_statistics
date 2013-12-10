#!/usr/bin/env Rscript

# generate_summary_statistics.R
# Generate dXY and partitioning statistics for models generated using run_model_combinations.pl

# Written for "Evaluating statistics for the identification of introgressed loci"
# by Simon H. Martin, John W. Davey and Chris D. Jiggins
# John Davey:   jd626@cam.ac.uk
# Simon Martin: shm45@cam.ac.uk
# November-December 2013

library(parallel)
suppressMessages(library(reshape))
library(plyr)
library(optparse)
library(stringr)
library(tools)

options<-list(make_option(c("-T","--threads"),default=1,help="Number of threads [default %default]",metavar="numthreads"),
              make_option(c("-m","--modelfiles"),default="",type="character",help="Folder of CSV files with model output"),
              make_option(c("-r","--realdata"),default="",type="character",help="CSV file containing windows from real genome"),
              make_option(c("-l","--listofmodels"),default="",type="character",help="CSV file of model parameter combinations")
         )

opt<-parse_args(OptionParser(option_list=options))

if (opt$listofmodels == "") stop("Please specify a list of model parameter combinations in CSV format with -l")
if (opt$modelfiles == "" & opt$realdata == "") stop("Please specify a folder of model files with -m or a CSV file containing real data with -r")
if (opt$modelfiles != "" & !file.exists(opt$modelfiles)) stop(paste("Directory", opt$modelfiles, "does not exist! Please specify a valid directory with -m"))
if (opt$realdata != "" & !file.exists(opt$realdata)) stop(paste("Real data file", opt$realdata, "does not exist! Please specify a valid file with -r"))

partition.x<-function(model,stat,threshold=0.1,win) {
    outnum<-threshold*win
    Partition<-sapply(rank(-model[,stat]), function(x) if (x<=outnum) "_Outlier" else "_Background")
    Partition<-paste0(stat,Partition)
    model<-cbind(model,Partition)
    names(model)[names(model)=="Partition"]<-paste0(stat,"_Partition")
    model
}

run.tests<-function(values,groups,test.name,test.func) {
    # Assumes order of groups is Background, Outlier
    formula<-values~groups
    two.tailed<-test.func(formula,alternative="two.sided",conf.int=TRUE)
    if (length(two.tailed$estimate)==1) mean.diff<-two.tailed$estimate else mean.diff<-two.tailed$estimate[1] - two.tailed$estimate[2]
    out<-data.frame(
        mean.diff=mean.diff,
        p.two.tailed=two.tailed$p.value,
        p.background.higher=test.func(formula,alternative="greater")$p.value,
        p.background.lower=test.func(formula,alternative="less")$p.value
    )
    names(out)<-paste(test.name,names(out),sep=".")
    out
}

summarise.tests<-function(model.sum,model.stats,model.partition,test.name,test.func) {
    model.tests<-t(apply(model.stats,2,function(x) unlist(run.tests(x,model.partition,test.name,test.func))))
    model.tests<-cbind(Variable=rownames(model.tests),model.tests)
    merge(model.sum,model.tests,by="Variable",sort=FALSE)
}

get.model.summary<-function(models,column,pars=NULL) {
    model.stats<-models[c("P1P2_dxy","P1P3_dxy","P2P3_dxy")]
    model.partition<-models[,column]
    
    model.means<-melt(aggregate(model.stats,list(model.partition),mean,na.rm=TRUE),id="Group.1")
    model.sds<-melt(aggregate(model.stats,list(model.partition),sd,na.rm=TRUE),id="Group.1")
    model.ses<-melt(aggregate(model.stats,list(model.partition), function(x) {x<-na.omit(x);sqrt(var(x)/length(x))}),id="Group.1")
    model.sum<-cbind(pars,model.means,model.sds$value,model.ses$value)
    names(model.sum)<-c(names(pars),"Model","Variable","Mean","SD","SE")

    model.sum<-summarise.tests(model.sum,model.stats,model.partition,"t",t.test)
    model.sum<-summarise.tests(model.sum,model.stats,model.partition,"w",wilcox.test)
    
    model.sum
}

calc.mfD0<-function(df) {
    apply(df,1,function(x) if ((is.na(x["D"])) || (as.numeric(x["D"])<0)) NA else as.numeric(x["mf"]))
}


summarise.model<-function(filename) {
    pars<-data.frame(str_match(filename,"Alternate_t123-(.+)_t23-(.+)\\.Background_t123-(.+)_t21-(.+)\\.(.+)\\.csv"),stringsAsFactors=FALSE)
    names(pars)<-c("File","Alternate_t123","Alternate_t23","Background_t123","Background_t21","Windows")
    pars$Windows<-as.numeric(pars$Windows)
    
    models<-read.csv(filename)

    models$Model<-factor(str_match(models$Model,"(.+)_t123")[,2],levels=c("Background","Alternate"))
    models$Model<-revalue(models$Model,c(Background="Real_Background",Alternate="Real_Outlier"))

    models$mfD0<-calc.mfD0(models)

    models<-partition.x(models,"D",0.1,pars$Windows)
    models<-partition.x(models,"mfD0",0.1,pars$Windows)

    model.summary<-get.model.summary(models,"Model",pars)
    D.summary<-get.model.summary(models,"D_Partition",pars)
    mfD0.summary<-get.model.summary(models,"mfD0_Partition",pars)
    
    model.sum<-rbind(model.summary,D.summary,mfD0.summary)
    
    is.D.positive<-models$D>=0
    model.comp<-data.frame(pars,
        Stat=c("D","mfD0"),
        Stat.Mean=c(mean(models$D,na.rm=TRUE)*100,mean(models$mfD0,na.rm=TRUE)*100),
        Stat.SD=c(sd(models$D,na.rm=TRUE)*100,sd(models$mfD0,na.rm=TRUE)*100),
        Dpos=sum(is.D.positive,na.rm=TRUE),
        DposAlternate=sum(is.D.positive & models$Model=="Real_Outlier",na.rm=TRUE),
        DposOutlier=c(sum(is.D.positive & models$D_Partition=="D_Outlier", na.rm=TRUE),
                      sum(is.D.positive & models$mfD0_Partition=="mfD0_Outlier", na.rm=TRUE)),
        OutlierAlternate=c(sum(models$Model=="Real_Outlier" & models$D_Partition=="D_Outlier",na.rm=TRUE),
                            sum(models$Model=="Real_Outlier" & models$mfD0_Partition=="mfD0_Outlier",na.rm=TRUE))
    )
    model.comp$OutlierAlternatePC<-model.comp$OutlierAlternate/model.comp$DposOutlier*100

    list(model.sum,model.comp)
}

options(width=200)

if (opt$modelfiles != "") {
    csv.summary<-mclapply(dir(opt$modelfiles,"*.csv",full.names=TRUE),summarise.model,mc.cores=opt$threads)

    models<-rbind.fill(lapply(csv.summary,function(x) x[[1]]))
    model.comp<-rbind.fill(lapply(csv.summary,function(x) x[[2]]))

    model.list<-read.csv(opt$listofmodels,stringsAsFactors=FALSE)

    models<-merge(models,unique(model.list[,c("Background_t123","Background_t21","Alternate_t123","Alternate_t23","ModelType")]),by=c("Background_t123","Background_t21","Alternate_t123","Alternate_t23"))
    write.table(models,file=paste0(opt$modelfiles,".dxy.summary.tsv"),sep="\t",row.names=FALSE)

    model.comp<-merge(model.comp,unique(model.list[,c("Background_t123","Background_t21","Alternate_t123","Alternate_t23","ModelType")]),by=c("Background_t123","Background_t21","Alternate_t123","Alternate_t23"))
    write.table(model.comp,file=paste0(opt$modelfiles,".partition.summary.tsv"),sep="\t",row.names=FALSE)
}

if (opt$realdata != "") {
    real<-read.csv(opt$realdata)
    real$mfD0<-calc.mfD0(real)
    real<-partition.x(real,"D",0.1,nrow(real))
    real<-partition.x(real,"mfD0",0.1,nrow(real))
    
    pars<-data.frame(File=opt$realdata,Windows=nrow(real))
    real.stats<-rbind(get.model.summary(real,"D_Partition",pars),get.model.summary(real,"mfD0_Partition",pars))
    write.table(real.stats,file=paste0(file_path_sans_ext(basename(opt$realdata)),".dxy.summary.tsv"),sep="\t",row.names=FALSE)
}