#!/usr/bin/env Rscript

# Figure_S2.R
# Simulation to examine the effect of recombination rate on the variance of D and f estimators - Figure S2

# Written for "Evaluating statistics for the identification of introgressed loci"
# by Simon H. Martin, John W. Davey and Chris D. Jiggins
# John Davey:   jd626@cam.ac.uk
# Simon Martin: shm45@cam.ac.uk
# August 2014



library(phyclust)
library(parallel)
library(plyr)



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



seqgen.to.DNAbin<-function(simseqs, len ) {

    bases<-matrix(nrow=length(simseqs),ncol=len,byrow=TRUE)
    rownames(bases)<-1:length(simseqs)
    for (s in 1:length(simseqs)) {
        seqstats<-unlist(strsplit(simseqs[s],"\\s+"))
        seqstats[1]<-gsub('s','',seqstats[1])
        bases[as.numeric(seqstats[1]),]<-unlist(strsplit(tolower(seqstats[2]),''))
    }
    as.DNAbin(bases)
}


get.pop.slice<-function(popsites,bi.sites) {
    counts<-apply(popsites,2,base.freq,freq=TRUE)
    alleles<-apply(counts/counts,2,sum,na.rm=TRUE)
    s<-length(which(alleles>1))
    list(counts=counts,counts.bi=counts[,bi.sites],alleles=alleles,s=s)
}



get.abba.baba.stats<-function(counts,bi,P1,P2,P3,P4,P31=NULL,P32=NULL) {
    
    do.real.f <- (is.null(P31)==F & is.null(P32)==F) #if both P3 subsets are given, calculate real f
    
    ABBA = 0
    BABA = 0
    maxABBA_G = 0
    maxBABA_G = 0
    maxABBA_hom = 0
    maxBABA_hom = 0
    maxABBA_D = 0
    maxBABA_D = 0

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
      if (do.real.f == T){
        P31df <- P31$counts.bi[derived,i] / sum(P31$counts.bi[,i])
        P32df <- P32$counts.bi[derived,i] / sum(P32$counts.bi[,i])
        }
      ABBA <- sum(ABBA, (1 - P1df) * P2df * P3df * (1 - P4df), na.rm = TRUE)
      BABA <- sum(BABA, P1df * (1 - P2df) * P3df * (1 - P4df), na.rm = TRUE)
      if (do.real.f==T){
        maxABBA_G <- sum(maxABBA_G, (1 - P1df) * P31df * P32df * (1 - P4df), na.rm = TRUE)
        maxBABA_G <- sum(maxBABA_G, P1df * (1 - P31df) * P32df * (1 - P4df), na.rm = TRUE)
          }
      maxABBA_hom <- sum(maxABBA_hom, (1 - P1df) * P3df * P3df * (1 - P4df), na.rm = TRUE)
      maxBABA_hom <- sum(maxBABA_hom, P1df * (1 - P3df) * P3df * (1 - P4df), na.rm = TRUE)
      if (is.na(P3df) == FALSE & is.na(P2df) == FALSE & P3df >= P2df){
        maxABBA_D <- sum(maxABBA_D, (1 - P1df) * P3df * P3df * (1 - P4df), na.rm = TRUE)
        maxBABA_D <- sum(maxBABA_D, P1df * (1 - P3df) * P3df * (1 - P4df), na.rm = TRUE)
        } else{
        maxABBA_D <- sum(maxABBA_D, (1 - P1df) * P2df * P2df * (1 - P4df), na.rm = TRUE)
        maxBABA_D <- sum(maxBABA_D, P1df * (1 - P2df) * P2df * (1 - P4df), na.rm = TRUE)
      }
    }
    output<-data.frame(ABBA=ABBA,BABA=BABA,
               maxABBA_G=maxABBA_G, maxBABA_G=maxBABA_G,
               maxABBA_D=maxABBA_D, maxBABA_D=maxBABA_D,
               D=(ABBA - BABA) / (ABBA + BABA),
               fG=(ABBA - BABA) / (maxABBA_G - maxBABA_G),
               fhom=(ABBA - BABA) / (maxABBA_hom - maxBABA_hom),
               fd=(ABBA - BABA) / (maxABBA_D - maxBABA_D)
               )
    output$fdD0<-apply(output,1,function(x) if ((is.na(x["D"])) || (as.numeric(x["D"]) < 0)) NA else as.numeric(x["fd"]))
    output
}



# Requires DNAbin input, vector of population names for each sequence in the sample,

# The arguments 'P1', 'P2', 'P3' and 'P4' (and optional pops 'P31' and 'P32') are
# each given as an integer vector of all sequences in each pop.
# eg. P1 = c(1,2,3,4,5,6,7,8), P2 = 9:16, P3 = 17:24, P4 = 25:32, P31 = 17:20, P32 = 21:24

ABBABABA <- function(aln, P1, P2, P3, P4, P31 = NULL, P32 = NULL) {
  

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

  P1.slice<-get.pop.slice(aln.var[P1,],bi.var)
  P2.slice<-get.pop.slice(aln.var[P2,],bi.var)
  P3.slice<-get.pop.slice(aln.var[P3,],bi.var)
  P4.slice<-get.pop.slice(aln.var[P4,],bi.var)
  if (is.null(P31)==F & is.null(P32)==F) {
    P31.slice<-get.pop.slice(aln.var[P31,],bi.var)
    P32.slice<-get.pop.slice(aln.var[P32,],bi.var)
    }

  abba.baba.stats<-get.abba.baba.stats(counts.all,bi.all,P1.slice,P2.slice,P3.slice,P4.slice,P31.slice,P32.slice)

  data.frame(S=s.all, P1_S=P1.slice$s, P2_S=P2.slice$s, P3_S=P3.slice$s, P4_S=P4.slice$s, abba.baba.stats)
}


sim_D_f <- function(n=1, nsam, P1, P2, P3, P4, P31, P32, ms_opts, len, scale, subModel, stats){
  ms.out<-ms(nsam=nsam, opts=ms_opts)
  ms.bases<-ms.to.DNAbin(ms.out, nsam)
  
  ms.trees<-ms.out[3:(grep("^segsites",ms.out)-1)]
  
  #simseqs<-system(sprintf("seq-gen -m%s -l %d -s %f -p %d -q", "HKY", len, scale, length(ms.trees)), input=ms.trees, intern=TRUE)[-1]
  
  seqgen.tmp.file<-paste("tmp", "seqgen", as.character(runif(1)), "out", sep=".")
  system(sprintf("seq-gen -m%s -l %d -s %f -p %d -q > %s", subModel, len, scale, length(ms.trees),  seqgen.tmp.file), input=ms.trees, ignore.stderr=TRUE)
  simseqs<-scan(seqgen.tmp.file, what="", sep="\n", quiet=TRUE)[-1]
  file.remove(seqgen.tmp.file)
  
  sg.bases<-seqgen.to.DNAbin(simseqs, len)
  
  abdf <- rbind(ABBABABA(ms.bases, P1,P2,P3,P4,P31,P32)[stats], ABBABABA(sg.bases, P1,P2,P3,P4,P31,P32)[stats])
  
  return(abdf[stats])
  }



#########################################################################################################################
### simulation. Multiple fs, multiple times of gene flow, using both ms and seqgen


subModel = "HKY"

P1 <- 1:8
P2 <- 9:16
P3 <- 17:24
P4 <- 25:32
P31 <- 17:20
P32 <- 21:24

nsam = 32

f = 0.1

nreps = 100

stats <- c("D","fG", "fhom","fd")

t = 50
len = 5000
rec_number = c(5,50,500)
tGF = 0.1
scale=0.01
theta = paste("-t", t)
pop_sizes = "-I 4 8 8 8 8"
join1 = "-ej 1 2 1"
join2 = "-ej 2 3 1"
join3 = "-ej 3 4 1"


### gene flow 3-2
# three dimensional array for the output. Dimensions are fs (0.1,0.2 etc), stats (D, f etc), and reps (1:100).
all_data_32 <- array(dim = c(nreps, length(rec_number), length(stats)), dimnames = list(1:nreps, c(0.001,0.01,0.1), stats))

#for each f we will generate nreps alignments, and get all the stats for each
for (x in 1:length(rec_number)){
  
  print(x)
  
  receiver = paste("-es",tGF, "2", 1 - f) #f is 1 - proportion staying. (See ms documentation)
  donor = paste("-ej", tGF, "5 3")
  
  rec = paste("-r", rec_number[x], len)
  
  ms_opts = paste("-T",theta,pop_sizes,join1,join2,join3,receiver,donor,rec)
  
  sim_data_list <- mclapply(1:nreps, sim_D_f, nsam=nsam, P1=P1, P2=P2, P3=P3, P4=P4, P31=P31, P32=P32, ms_opts=ms_opts, len=len, scale=scale, subModel=subModel, stats=stats, mc.cores = 3)
  
  sim_data_array <- aaply(laply(sim_data_list, as.matrix),c(1,3),function(x) x) 
  
  all_data_32[,x,] <- sim_data_array[,,2]
          
  }
  

all_data_32 <- ifelse(abs(all_data_32 == "Inf"), NA, all_data_32) # remove all values estimated calculated to be Inf or -Inf



######################################################################################
### plot 

png(file = paste("Figure_S2.png", sep = ""), width = 4000, height = 1500, res = 400)

par(mfrow = c(1,4), mar = c(3.5,3,2,1))

boxplot(all_data_32[,,"D"], ylim = c(-1,1), main = expression(italic("D")))
mtext(side = 1, text = "recombination rate (4Nr)", line = 2.3)


boxplot(all_data_32[,,"fG"], ylim = c(-1,1), main = expression(italic("f"["G"])))
mtext(side = 1, text = "recombination rate (4Nr)", line = 2.3)


boxplot(all_data_32[,,"fhom"], ylim = c(-1,1), main = expression(italic("f"["hom"])))
mtext(side = 1, text = "recombination rate (4Nr)", line = 2.3)


boxplot(all_data_32[,,"fd"], ylim = c(-1,1), main = expression(italic("f"["d"])))
mtext(side = 1, text = "recombination rate (4Nr)", line = 2.3)


dev.off()


