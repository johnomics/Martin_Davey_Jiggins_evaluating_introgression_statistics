#!/usr/bin/env Rscript

# Figures_3_S2_S3.R
# Plots for whole genome Heliconius data

# Written for "Evaluating statistics for the identification of introgressed loci"
# by Simon H. Martin, John W. Davey and Chris D. Jiggins
# John Davey:   jd626@cam.ac.uk
# Simon Martin: shm45@cam.ac.uk
# May 2014


############## Functions ###############

#function to rearrange scaffolds by chromosome and add a chromosome column
add.chrom <- function(table,agp) {
  chromosome <- vector(length = length(table[,1]))
  for (x in 1:length(table[,1])){
    scaf <- table$scaffold[x]
    if (table$scaffold[x] %in% agp$scaffold){
      y <- which(agp$scaffold == scaf)
      chromosome[x] <- agp$chromosome[y]
      }
    else{
      chromosome[x] <- NA
      }
    }
  return(cbind(chromosome,table))
  }



########### input data ################

#agp file carrys chromosome information for mapped scaffolds
# Available in Data Dryad repository http://datadryad.org/resource/doi:10.5061/dryad.dk712
agp <- read.delim("Hmel1-1_hox_RAD_matepair_chromosomes_Zupdated.agp", as.is = T, sep = "\t", header = F)

names(agp) = c("chromosome", "start", "end", "order", "DN", "scaffold", "one", "length", "ori", "rec")
agp <- subset(agp, DN == "D")

#remove the word "_unmapped" from chromnames)
for (chr in c(1:20, "Z")) {
    agp$chromosome <- ifelse(agp$chromosome == paste0("chr",chr,"_unmapped"), paste0("chr",chr), agp$chromosome)
}

#input data


ABBA_table <- read.csv("Heliconius_autosome_windows.csv", as.is = T)
ABBA_table_Z <- read.csv("Heliconius_Zchromosome_windows.csv", as.is = T)
ABBA_table <- rbind(ABBA_table, ABBA_table_Z)

#filter for windows with enough sites above the minimum data cutoff
ABBA_table <- subset(ABBA_table, sitesOverMinExD >= 3000)


#calculate mean Pi in the three pops

mean_Pi <- apply(ABBA_table[,c("P1_pi", "P2_pi", "P3_pi")], 1, mean, na.rm = T) 

ABBA_table <- cbind(ABBA_table, mean_Pi)


#scaffolds containing wing patterning loci
ybsub <- subset(ABBA_table, scaffold == "HE667780" & position >= 650000 & position <= 900000)
bdsub <- subset(ABBA_table, scaffold == "HE670865" & position >= 300000 & position <= 450000)


###### plotting D and f against number of segregating sites


png(file = "Figure_3.png", width = 1500, height = 4500, res = 500)


par(mfrow = c(3,1), mar = c(3,3,1,1))

plot(mean_Pi, ABBA_table$D, cex = 1,  , pch = 21, col = rgb(0,0,0,0), bg = rgb(0,0,0,0.03), bty = "n", xlab = "", ylab = "", xlim = c(0,0.08) )
points(ybsub$mean_Pi, ybsub$D, cex = 1, pch = 1, col = rgb(1,0.8,0,1))
points(bdsub$mean_Pi, bdsub$D, cex = 1, pch = 1, col = "red")


plot(ABBA_table$mean_Pi[ABBA_table$D >= 0], ABBA_table$fd[ABBA_table$D >= 0], cex = 1,  , pch = 21, col = rgb(0,0,0,0), bg = rgb(0,0,0,0.03), ylim = c(0,1), bty = "n", xlab = "", ylab = "", xlim = c(0,0.08) )
points(ybsub$mean_Pi[ybsub$D >= 0], ybsub$fd[ybsub$D >= 0], cex = 1, pch = 1, col = rgb(1,0.8,0,1))
points(bdsub$mean_Pi[bdsub$D >= 0], bdsub$fd[bdsub$D >= 0], cex = 1, pch = 1, col = "red")



### variance by number of segregating sites - plotted by chromosome
library(RColorBrewer)
cols <- c(brewer.pal(12,"Set3"),brewer.pal(9,"Set1"))

#rearrange scaffolds by chromosome (this is only necessary for the third, chromosomal, plot)
WG_table <- add.chrom(ABBA_table,agp)


chromNames <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chrZ")


#means and variances
vD <- numeric(length <- length(chromNames))
vfd <- numeric(length <- length(chromNames))
mean_mean_Pi <- numeric(length <- length(chromNames))

for (x in 1:length(chromNames)){
  vD[x] <- var(WG_table$D[WG_table$chromosome == chromNames[x]], na.rm = T)
  vfd[x] <- var(WG_table$fd[WG_table$chromosome == chromNames[x] & WG_table$D >= 0], na.rm = T)
  mean_mean_Pi[x] <- mean(WG_table$mean_Pi[WG_table$chromosome == chromNames[x]], na.rm = T)
  }


plot(vD ~ mean_mean_Pi, ylim = c(0,0.16), xlim = c(0.01,0.025), pch = 19, cex = 1.5, col = cols, xlab = "", ylab = "", bty = "n")

plot(vfd ~ mean_mean_Pi, pch = 15, cex = 1.5, col = cols)

legend(0.02, 0.15,legend = c(1:20,"Z"), pch = 19, col = cols, ncol = 3, bty = "n", cex = 0.8, pt.cex = 1.1)

dev.off()


#### D and all fs agains pi


png(file = "Figure_S2.png", width = 3500, height = 3500, res = 400)


par(mfrow = c(2,2), mar = c(3.5,3.5,1,1))

plot(ABBA_table$mean_Pi, ABBA_table$D, cex = 1,  , pch = 21, col = rgb(0,0,0,0), bg = rgb(0,0,0,0.03), bty = "n", xlab = "", ylab = "", xlim = c(0,0.08) )
mtext(side=1, text = "Nucleotide diversity", line = 2)
mtext(side=2, text = expression(italic("D")), line = 2)
points(ybsub$mean_Pi, ybsub$D, cex = 1, pch = 1, col = rgb(1,0.8,0,1))
points(bdsub$mean_Pi, bdsub$D, cex = 1, pch = 1, col = "red")


plot(ABBA_table$mean_Pi[ABBA_table$D >= 0], ABBA_table$fG[ABBA_table$D >= 0], cex = 1,  , pch = 21, col = rgb(0,0,0,0), bg = rgb(0,0,0,0.03), ylim = c(0,1), bty = "n", xlab = "", ylab = "", xlim = c(0,0.08) )
mtext(side=1, text = "Nucleotide diversity", line = 2)
mtext(side=2, text = expression(italic("f"["G"])), line = 2)
points(ybsub$mean_Pi[ybsub$D >= 0], ybsub$fG[ybsub$D >= 0], cex = 1, pch = 1, col = rgb(1,0.8,0,1))
points(bdsub$mean_Pi[bdsub$D >= 0], bdsub$fG[bdsub$D >= 0], cex = 1, pch = 1, col = "red")

plot(ABBA_table$mean_Pi[ABBA_table$D >= 0], ABBA_table$fhom[ABBA_table$D >= 0], cex = 1,  , pch = 21, col = rgb(0,0,0,0), bg = rgb(0,0,0,0.03), ylim = c(0,1), bty = "n", xlab = "", ylab = "", xlim = c(0,0.08) )
mtext(side=1, text = "Nucleotide diversity", line = 2)
mtext(side=2, text = expression(italic("f"["hom"])), line = 2)
points(ybsub$mean_Pi[ybsub$D >= 0], ybsub$fhom[ybsub$D >= 0], cex = 1, pch = 1, col = rgb(1,0.8,0,1))
points(bdsub$mean_Pi[bdsub$D >= 0], bdsub$fhom[bdsub$D >= 0], cex = 1, pch = 1, col = "red")

plot(ABBA_table$mean_Pi[ABBA_table$D >= 0], ABBA_table$fd[ABBA_table$D >= 0], cex = 1,  , pch = 21, col = rgb(0,0,0,0), bg = rgb(0,0,0,0.03), ylim = c(0,1), bty = "n", xlab = "", ylab = "", xlim = c(0,0.08) )
mtext(side=1, text = "Nucleotide diversity", line = 2)
mtext(side=2, text = expression(italic("f"["d"])), line = 2)
points(ybsub$mean_Pi[ybsub$D >= 0], ybsub$fd[ybsub$D >= 0], cex = 1, pch = 1, col = rgb(1,0.8,0,1))
points(bdsub$mean_Pi[bdsub$D >= 0], bdsub$fd[bdsub$D >= 0], cex = 1, pch = 1, col = "red")


dev.off()








### variances in mean_Pi bins

pi_cuts <- cut(mean_Pi, seq(0,0.06,0.01))

D_vars <- sapply(1:length(levels(pi_cuts)), function(x){var(ABBA_table$D[pi_cuts == levels(pi_cuts)[x] & ABBA_table$D >= 0], na.rm = T)})
fd_vars <- sapply(1:length(levels(pi_cuts)), function(x){var(ABBA_table$fd[pi_cuts == levels(pi_cuts)[x] & ABBA_table$D >= 0] , na.rm = T)})
fG_vars <- sapply(1:length(levels(pi_cuts)), function(x){var(ABBA_table$fG[pi_cuts == levels(pi_cuts)[x]  & ABBA_table$D >= 0 & ABBA_table$fG >= 0 & ABBA_table$fG <= 1], na.rm = T)})
fhom_vars <- sapply(1:length(levels(pi_cuts)), function(x){var(ABBA_table$fhom[pi_cuts == levels(pi_cuts)[x] & ABBA_table$D >= 0 & ABBA_table$fhom >= 0 & ABBA_table$fhom <= 1], na.rm = T)})

png("Figure_S3.png", width = 2500, height = 2500, res = 400)
par(mar = c(4,4,1,1))

plot(D_vars, type = "b", pch = 19, xlab = "Pi bin", xaxt = "n", ylab = "Variance", ylim = c(0,0.1))
axis(1, at=1:6, labels = c("0 - 0.01", "0.01 - 0.02", "0.02 - 0.03", "0.03 - 0.04", "0.04 - 0.05", "0.05 - 0.06"), las = 1) 
points(fd_vars, type = "b")
points(fG_vars, type = "b", pch = 2)
points(fhom_vars, type = "b", pch = 3)

legend(5,0.1, legend = c(expression(italic("D")), expression(italic("f"["G"])), expression(italic("f"["hom"])), expression(italic("f"["d"]))), pch = c(19,1,2,3), bty = "n")

dev.off()

