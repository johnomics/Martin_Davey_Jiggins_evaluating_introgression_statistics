#!/usr/bin/env Rscript

# shared_ancestry_simulator.R
# Generate models using ms and seq-gen based on model specified by YAML file
# Requires seq-gen to be installed separately

# Written for "Evaluating the use of ABBA-BABA statistics to locate introgressed loci"
# by Simon H. Martin, John W. Davey and Chris D. Jiggins
# Simon Martin: shm45@cam.ac.uk
# John Davey:   jd626@cam.ac.uk
# October-December 2013, May 2014, August 2014

library(tools)
library(parallel)
library(plyr)
library(yaml)
library(optparse)
library(ggplot2,quietly=TRUE)
library(gridExtra,quietly=TRUE)
suppressMessages(library(reshape,quietly=TRUE))

# phyclust contains ms,
# but also depends on ape,
# making the DNAbin class available
library(phyclust, quietly=TRUE)

isnt.null <- function(x)!is.null(x)
isnt.na <- function(x)!is.na(x)
se <- function(x) {sd(x, na.rm = T)/sqrt(length(which(isnt.na(x))))}

# Parse options
options<-list(make_option(c("-t","--threads"),default=1,help="Number of threads [default %default]",metavar="numthreads"),
              make_option(c("-w","--windows"),default=100,help="Number of windows to simulate",metavar="numwindows"),
              make_option(c("-c","--combined"),type="character",help="YAML file list with proportions for combined models, eg -c model1.yml:0.1,model2.yml:0.9"),
              make_option(c("-v","--verbose"),action="store_true", default=FALSE, help="Print ms command to standard out")
         )

opt<-parse_args(OptionParser(option_list=options))

if (is.null(opt$combined)) stop("Please specify model files in YAML format and proportions with -c, eg -c model1.yml:0.1,model2.yml:0.9")
if(!is.numeric(opt$windows) | opt$windows<1) stop("Please specify a positive number of windows with -w")
if(!is.numeric(opt$threads) | opt$threads<1) stop("Please specify a positive number of threads with -T")

ms.to.DNAbin <- function(ms_output, nsam=0, inc_inv=FALSE, len=NULL) {
    snp.strings <- ms_output[(length(ms_output)-nsam+1):length(ms_output)]
    binary.mat.var <- t(matrix(as.numeric(unlist(sapply(snp.strings,strsplit, split = ""))),ncol = length(snp.strings)))
    if (inc_inv==TRUE) {
        stopifnot(mode(len)=="numeric" & len >= length(binary.mat.var[1,]))
        nsam <- length(binary.mat.var[,1])
        nvar <- length(binary.mat.var[1,])
        binary.mat <- matrix(numeric(length=nsam*len), nrow=nsam)
        positions <- floor(as.numeric(unlist(strsplit(ms_output[length(ms_output)-nsam],split="    "))[2:(nvar+1)])*len)+1
        binary.mat[1:nsam,positions] <- binary.mat.var
    } else {
        binary.mat <- binary.mat.var
    }
    bases <- ifelse(binary.mat == "1", "t", "a")
    rownames(bases)<-1:nsam
    as.DNAbin(bases)
}


#Requires DNAbin input and vector of population names for each sequence in the sample
nucdiv <- function(aln, popIndex) {
  #check if pops were defined
  if (is.null(popIndex) == TRUE) {
    stop("")
    pops <- list( "1" = labels(aln))
    } else {
    pops <- list()
    for (popName in unique(popIndex)) {
      pops[[popName]] <- labels(aln)[which(popIndex == popName)]
    }
  }
  #create an output matrix. It will have all populations as rows and columns, so upper and lower triangles will be mirrored.
  #this just makes indexing easier. Each block will be dxy for that pair pop populations, except diagnols which wll be pi.
  output <- matrix(nrow = length(pops), ncol = length (pops))
  rownames(output) <- names(pops)
  colnames(output) <- names(pops)
  #Get pairwise distance matrix 
  pdm <- dist.dna(aln, model = "raw", pairwise.deletion=TRUE, as.matrix = TRUE)
  #Now, for each pop, get pi and for each pair get dxy
  #Pi is the mean of the lower triangle, exclusing the diagnol, when the same pop is used to index the rows and columns:
  for (X in 1:length(pops)){
    for (Y in X:length(pops)){
      sub_pdm <- pdm[pops[[X]],pops[[Y]]]
      #remove upper triangle and, if they're the same pop, remove diagnol, then take the mean and add to output matrix
      if (X == Y){
        sub_pdm[upper.tri(sub_pdm, diag = TRUE)] <- NA
        output[X,Y] <- mean(sub_pdm, na.rm = T)
      } else {
        sub_pdm[upper.tri(sub_pdm, diag = FALSE)] <- NA
        output[X,Y] <- mean(sub_pdm, na.rm = T)
        output[Y,X] <- mean(sub_pdm, na.rm = T) # add to both lower and upper - for ease of indexing
      }
    }
  }
  output
}


nucdiv4pop<-function(aln,popIndex) {
    if (length(unique(popIndex)) != 4) stop("Please supply four populations only for nucleotide divergence calculations")
    ndmat<-nucdiv(aln,popIndex)
    data.frame(
        P1_Pi=ndmat[1,1],P2_Pi=ndmat[2,2],P3_Pi=ndmat[3,3],P4_Pi=ndmat[4,4],
        P1P2_dxy=ndmat[1,2],P1P3_dxy=ndmat[1,3],P1P4_dxy=ndmat[1,4],P2P3_dxy=ndmat[2,3],P2P4_dxy=ndmat[2,4],P3P4_dxy=ndmat[3,4]
    )
}


get.pop.slice<-function(popsites,bi.sites) {
    counts<-apply(popsites,2,base.freq,freq=TRUE)
    alleles<-apply(counts/counts,2,sum,na.rm=TRUE)
    s<-length(which(alleles>1))
    list(counts=counts,counts.bi=counts[,bi.sites],alleles=alleles,s=s)
}


get.abba.baba.stats<-function(counts,bi,P1,P2,P3,P4,P31=NULL,P32=NULL) {

    ABBA <- 0
    BABA <- 0
    maxABBA_G <- 0
    maxBABA_G <- 0
    maxABBA_hom <- 0
    maxBABA_hom <- 0
    maxABBA_D <- 0
    maxBABA_D <- 0

    do.real.f <- (isnt.null(P31) & isnt.null(P32)) #if both P3 subsets are given, calculate real f

    for (i in 1:length(bi)){
        i.bi.counts<-counts[,bi[i]]
        alleles <- names(which(i.bi.counts > 0))
        anc <- names(which(P4$counts.bi[,i] != 0 ))
        if (length(anc) != 1) { anc <- names(which(i.bi.counts == max(i.bi.counts)))[1] }
        derived = alleles[which(alleles != anc)]

        P1df <- P1$counts.bi[derived,i] / sum(P1$counts.bi[,i])
        P2df <- P2$counts.bi[derived,i] / sum(P2$counts.bi[,i])
        P3df <- P3$counts.bi[derived,i] / sum(P3$counts.bi[,i])
        P4df <- P4$counts.bi[derived,i] / sum(P4$counts.bi[,i])
        
        if (do.real.f) {
          P31df <- P31$counts.bi[derived,i] / sum(P31$counts.bi[,i])
          P32df <- P32$counts.bi[derived,i] / sum(P32$counts.bi[,i])
        }
        
        ABBA <- sum(ABBA, (1 - P1df) * P2df * P3df * (1 - P4df), na.rm = TRUE)
        BABA <- sum(BABA, P1df * (1 - P2df) * P3df * (1 - P4df), na.rm = TRUE)
        
        if (do.real.f) {
            maxABBA_G <- sum(maxABBA_G, (1 - P1df) * P31df * P32df * (1 - P4df), na.rm = TRUE)
            maxBABA_G <- sum(maxBABA_G, P1df * (1 - P31df) * P32df * (1 - P4df), na.rm = TRUE)
        }
        
        maxABBA_hom <- sum(maxABBA_hom, (1 - P1df) * P3df * P3df * (1 - P4df), na.rm = TRUE)
        maxBABA_hom <- sum(maxBABA_hom, P1df * (1 - P3df) * P3df * (1 - P4df), na.rm = TRUE)
        
        if (isnt.na(P3df) & isnt.na(P2df) & P3df >= P2df) {
            maxABBA_D <- sum(maxABBA_D, (1 - P1df) * P3df * P3df * (1 - P4df), na.rm = TRUE)
            maxBABA_D <- sum(maxBABA_D, P1df * (1 - P3df) * P3df * (1 - P4df), na.rm = TRUE)
        } else {
            maxABBA_D <- sum(maxABBA_D, (1 - P1df) * P2df * P2df * (1 - P4df), na.rm = TRUE)
            maxBABA_D <- sum(maxBABA_D, P1df * (1 - P2df) * P2df * (1 - P4df), na.rm = TRUE)
        }
    }
    output<-data.frame(ABBA=ABBA,BABA=BABA,
               maxABBA_G=maxABBA_G, maxBABA_G=maxBABA_G,
               maxABBA_hom=maxABBA_hom, maxBABA_hom=maxBABA_hom,
               maxABBA_D=maxABBA_D, maxBABA_D=maxBABA_D,
               D=(ABBA - BABA) / (ABBA + BABA),
               fG=(ABBA - BABA) / (maxABBA_G - maxBABA_G),
               fhom=(ABBA - BABA) / (maxABBA_hom - maxBABA_hom),
               fd=(ABBA - BABA) / (maxABBA_D - maxBABA_D)
               )
    output$fdD0<-apply(output,1,function(x) if ((is.na(x["D"])) || (as.numeric(x["D"]) < 0)) NA else as.numeric(x["fd"]))
    output
}


ABBABABA <- function(aln, poplist) {
    #get allele counts to find variable sites
    counts.all<-apply(aln,2,base.freq,freq=TRUE)
    #get numbers of alleles at each site (counts/counts converts allele counts to 0 or 1)
    num.alleles.all<-colSums(counts.all/counts.all,na.rm=TRUE)
      
    #biallelic sites
    bi.all <- which(num.alleles.all == 2)

    #varying sites
    var.sites.all<-num.alleles.all>1
    aln.var<-aln[,var.sites.all]
    alleles.var<-num.alleles.all[var.sites.all]
    s.all<-length(alleles.var)
      
    bi.var<-which(alleles.var==2)


    P1.slice<-get.pop.slice(aln.var[poplist$P1,], bi.var)
    P2.slice<-get.pop.slice(aln.var[poplist$P2,], bi.var)
    P3.slice<-get.pop.slice(aln.var[poplist$P3,], bi.var)
    P4.slice<-get.pop.slice(aln.var[poplist$P4,], bi.var)

    P3.split.index<-floor(length(poplist$P3)/2)
    P31<-poplist$P3[1:P3.split.index]
    P32<-poplist$P3[(P3.split.index+1):length(poplist$P3)]
    P31.slice<-get.pop.slice(aln.var[P31,],bi.var)
    P32.slice<-get.pop.slice(aln.var[P32,],bi.var)

    abba.baba.stats<-get.abba.baba.stats(counts.all, bi.all,
                                         P1.slice,   P2.slice,
                                         P3.slice,   P4.slice,
                                         P31.slice,  P32.slice)
      
    data.frame(S=s.all, P1_S=P1.slice$s, P2_S=P2.slice$s, P3_S=P3.slice$s, P4_S=P4.slice$s, abba.baba.stats)
}

check.pops<-function(pops,totalpops,popnum) {
    if (length(pops) != popnum) stop(paste("Expected",popnum,"pops, got ",length(pops)))
    for (i in pops) {
        if (class(i) != "integer") stop("All population values must be integers")
        if (i > totalpops) stop("Population value larger than number of populations")
    }
}

parse.pops<-function(model) {
    if (!("pops" %in% names(model$ms))) stop("Please specify ms populations using the name 'pops'")

    model$ms$popnum<-length(model$ms$pops)
    model$ms$nsam<-0
    model$ms$popnames<-vector()
    model$ms$poplist<-list()
    model$ms$popsizes<-rep(NA,model$ms$popnum)
    for (pop in 1:model$ms$popnum) {
        model$ms$poplist[[pop]]<-(model$ms$nsam+1):(model$ms$nsam+model$ms$pops[[pop]]$seqs)
        model$ms$nsam<-model$ms$nsam + model$ms$pops[[pop]]$seqs
        model$ms$popnames<-c(model$ms$popnames,rep(model$ms$pops[[pop]]$name,model$ms$pops[[pop]]$seqs))
        model$ms$popsizes[pop]<-model$ms$pops[[pop]]$seqs
    }
    names(model$ms$poplist)<-paste0("P",1:model$ms$popnum)
    
    model
}

parse.migrations<-function(model) {

    if ("migration" %in% names(model$ms$opts)) {
        if (!("n0" %in% names(model$ms$opts))) stop("Found migration options but no N0; please specify n0 ms option")
        for (i in model$ms$opts$migration) {
            if (!("mij" %in% names(i))) stop("Please specify mij for all migration options")
            if (class(i$mij) != "numeric") stop("migration mij parameters should be numeric")

            if (!("pops" %in% names(i))) stop ("Please specify pops for all migration options")
            check.pops(i$pops,length(model$ms$pops),2)
        }
    }

    model
}

parse.splits<-function(model) {
    if ("split" %in% names(model$ms$opts)) {
        for (i in model$ms$opts$split) {
            if (!("time" %in% names(i)) & !("units" %in% names(i))) stop("Please specify time or units for all split options")
            if (!("pops" %in% names(i))) stop("Please specify pops for all split options")
            check.pops(i$pops,length(model$ms$pops),2)
        }
    }
    
    model
}

parse.n0<-function(model) {
    if ("n0" %in% names(model$ms$opts)) {
        n0test<-eval(parse(text=model$ms$opts$n0))
        if (class(n0test) != "numeric") stop("n0 must be a number or an expression producing a number")
        if (length(n0test) != 1) stop("n0 should only be a single number")
    }
    model
}

parse.theta<-function(model) {
    model$ms$opts$theta <- model$seqgen$s * model$seqgen$l
    model$ms$thetaopt<-paste("-t",model$ms$opts$theta)
    model
}

parse.recombination<-function(model) {
    if ("recombination" %in% names(model$ms$opts)) {
        if (!is.numeric(model$ms$opts$recombination)) stop("Please supply a number for rho (recombination)")
        model$ms$recombopt<-paste("-r",model$ms$opts$recombination,model$seqgen$l)
    } else {
        model$ms$recombopt <- ""
    }
    model
}

parse.flows<-function(model) {
    if ("flow" %in% names(model$ms$opts)) {
        for (i in model$ms$opts$flow) {
            if (!("time" %in% names(i))) stop("Please specify time for all flow options")
            if (!("f" %in% names(i))) stop("Please specify f (proportion of genome flowing) for all flow options")
            check.pops(i$pops, length(model$ms$pops),2)
        }
    }
    model
}

parse.seqgen<-function(model) {
    if (!("seqgen" %in% names(model))) stop("Please specify seqgen options using the name 'seqgen'")
    msg<-names(model$seqgen)
    if (!("s" %in% msg) && !("m" %in% msg) && !("l" %in% msg)) stop("Please specify options m, l and s for seqgen")
    return()
}

parse.model<-function(filename) {
    model<-yaml.load_file(filename)

    parse.seqgen(model)

    model<-parse.pops(model)

    if (!("opts" %in% names(model$ms))) stop("Please specify ms options using the name 'opts'")

    model<-parse.n0(model)
    model<-parse.theta(model)
    model<-parse.recombination(model)

    model<-parse.migrations(model)
    model<-parse.splits(model)
    model<-parse.flows(model)

    model$ms$optstring<-paste(c("-T","-I",model$ms$popnum,model$ms$popsizes),collapse=" ")
    if (model$ms$thetaopt != "") model$ms$optstring<-paste(model$ms$optstring,model$ms$thetaopt)
    if (model$ms$recombopt != "") model$ms$optstring<-paste(model$ms$optstring,model$ms$recombopt)
    model
}

get.migrations<-function(model) {
    if ("migration" %in% names(model$ms$opts)) {
        for (i in model$ms$opts$migration) {
            mig.rate<-i$mij * 4 * model$ms$opts$n0val
            migopt<-paste(c("-m",i$pops,mig.rate),collapse=" ")
            model$ms$optstring<-paste(model$ms$optstring,migopt)
        }
    }
    model
}

get.splits<-function(model) {
    if ("split" %in% names(model$ms$opts)) {
        for (i in model$ms$opts$split) {
            if ("units" %in% names(i)) {
                split.time<-i$units
            } else {
                split.time<-(i$time * 1000000 * model$ms$opts$genperyear) / (4 * model$ms$opts$n0val)
            }
            splitopt<-paste(c("-ej",split.time,i$pops),collapse=" ")
            model$ms$optstring<-paste(model$ms$optstring,splitopt)
        }
    }
    model
}

get.flows<-function(model) {
    if ("flow" %in% names(model$ms$opts)) {
        sparepop<-length(model$ms$pops)+1
        for (i in model$ms$opts$flow) {
            flow.optstring<-paste("-es",i$time,i$pops[2],i$f,"-ej",i$time,sparepop,i$pops[1])
            model$ms$optstring<-paste(model$ms$optstring,flow.optstring)
            sparepop<-sparepop+1
        }
    }

    model
}

get.model.options<-function(model) {

    model$ms$opts$n0val<-eval(parse(text=model$ms$opts$n0))

    model<-get.migrations(model)
    model<-get.splits(model)
    model<-get.flows(model)

    model
}

seqgen.to.DNAbin<-function(simseqs, model) {

    bases<-matrix(nrow=model$ms$nsam,ncol=model$seqgen$l,byrow=TRUE)
    rownames(bases)<-1:model$ms$nsam
    for (s in 1:model$ms$nsam) {
        seqstats<-unlist(strsplit(simseqs[s],"\\s+"))
        seqstats[1]<-gsub('s','',seqstats[1])
        bases[as.numeric(seqstats[1]),]<-unlist(strsplit(tolower(seqstats[2]),''))
    }
    as.DNAbin(bases)
}

get.seq.stats<-function(bases, prefix, model) {
    nddf<-nucdiv4pop(bases, model$ms$popnames)
    abdf<-ABBABABA(bases, model$ms$poplist)
    names(nddf)<-paste0(prefix, names(nddf))
    names(abdf)<-paste0(prefix, names(abdf))
    cbind(nddf, abdf)
}

simulate.window<-function(win, model) {
    model<-get.model.options(model)
    if (opt$verbose) print(model$ms$optstring)
    windf<-data.frame(blockNumber=win,blockName=win)
    pardf<-data.frame(scalerate=model$seqgen$s,length=model$seqgen$l)

    # Generate ms model
    ms.out<-ms(nsam=model$ms$nsam, opts=paste(model$ms$optstring, "-t", model$ms$opts$theta))
    # Get bases using ms variant sites
    ms.bases<-ms.to.DNAbin(ms.out, model$ms$nsam, inc_inv=TRUE, len=model$seqgen$l)

    # Get bases using seqgen
    ms.trees<-ms.out[3:(grep("^segsites",ms.out)-1)]
    seqgen.tmp.file<-paste("tmp", "seqgen", model$seqgen$m, model$seqgen$l, model$seqgen$s, win, "out", sep=".")
    system(sprintf("seq-gen -m%s -l %d -s %f -p %d -q > %s",
                            model$seqgen$m, model$seqgen$l, model$seqgen$s, length(ms.trees),
                            seqgen.tmp.file),
                            input=ms.trees,
                            ignore.stderr=TRUE)
    simseqs<-scan(seqgen.tmp.file, what="", sep="\n", quiet=TRUE)[-1]
    sg.bases<-seqgen.to.DNAbin(simseqs, model)
    file.remove(seqgen.tmp.file)
    
    cbind(windf,pardf,get.seq.stats(ms.bases, "ms.", model), get.seq.stats(sg.bases, "sg.", model))
}

simulate.model<-function(filename,proportion=1,cumprop=0) {
    propwin<-ceiling(opt$windows*proportion)
    startwin<-ceiling(opt$windows*cumprop)+1
    endwin<-startwin+propwin-1
    if (file.access(filename)==-1) stop(sprintf("Model file %s does not exist",filename))
    model<-parse.model(filename)
    rbind.fill(mclapply(startwin:endwin,function(win) simulate.window(win,model),mc.cores=opt$threads))
}

options(width=200)
propstats<-NULL
if (!(is.null(opt$combined))) {
    prop<-unlist(strsplit(opt$combined,","))
    modelfiles<-sapply(prop,function(p) unlist(strsplit(p,":"))[1])
    modelprops<-as.numeric(sapply(prop,function(p) unlist(strsplit(p,":"))[2]))
    opt$propdf<-data.frame(Filename=modelfiles,Proportion=modelprops,CumProp=cumsum(modelprops)-modelprops,stringsAsFactors=FALSE)
    if (sum(opt$propdf$Proportion) != 1) stop("Combined model proportions must add up to 1")

    propstats<-merge_all(mapply(simulate.model, opt$propdf$Filename,opt$propdf$Proportion,opt$propdf$CumProp,SIMPLIFY=FALSE))

    propmodel<-unlist(mapply(function(f,p) rep(f,p*opt$windows),opt$propdf$Filename,opt$propdf$Proportion))
    propstats<-cbind(propstats,Model=as.character(propmodel))

    ms.stats<-subset(propstats,select=grep("^sg.",names(propstats),invert=TRUE)) # Select all non-sg columns
    sg.stats<-subset(propstats,select=grep("^ms.",names(propstats),invert=TRUE)) # Select all non-ms columns
    names(ms.stats)<-gsub("^ms.","",names(ms.stats))
    names(sg.stats)<-gsub("^sg.","",names(sg.stats))

    propfilestub<-paste0(c(file_path_sans_ext(basename(opt$propdf$Filename)),opt$windows),collapse=".")
    write.csv(ms.stats,paste0(propfilestub,".ms.csv"))
    write.csv(sg.stats,paste0(propfilestub,".sg.csv"))
}
