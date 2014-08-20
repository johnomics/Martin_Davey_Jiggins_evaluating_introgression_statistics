#!/usr/bin/env Rscript

# Figure_S4.R
# Plots of variance for windows binned by pi, four different window sizes - Figure S4.

# Written for "Evaluating statistics for the identification of introgressed loci"
# by Simon H. Martin, John W. Davey and Chris D. Jiggins
# John Davey:   jd626@cam.ac.uk
# Simon Martin: shm45@cam.ac.uk
# August 2014



png("Figure_S4.png", width = 4000, height = 4000, res = 400)
par(mar = c(4,4,1,1), mfrow = c(2,2))


ABBA_table <- read.csv("Heliconius_autosome_windows_5kb.csv", as.is = T)
ABBA_table_Z <- read.csv("Heliconius_Zchromosome_windows_5kb.csv", as.is = T)
ABBA_table <- rbind(ABBA_table, ABBA_table_Z)

#filter for windows with enough sites above the minimum data cutoff
ABBA_table <- subset(ABBA_table, sitesOverMinExD >= 3000)


#calculate mean Pi in the three pops

mean_Pi <- apply(ABBA_table[,c("P1_pi", "P2_pi", "P3_pi")], 1, mean, na.rm = T) 

### variances in mean_Pi bins

pi_cuts <- cut(mean_Pi, seq(0,0.04,0.01))

D_vars <- sapply(1:length(levels(pi_cuts)), function(x){var(ABBA_table$D[pi_cuts == levels(pi_cuts)[x] & ABBA_table$D >= 0], na.rm = T)})
fd_vars <- sapply(1:length(levels(pi_cuts)), function(x){var(ABBA_table$fd[pi_cuts == levels(pi_cuts)[x] & ABBA_table$D >= 0] , na.rm = T)})
fG_vars <- sapply(1:length(levels(pi_cuts)), function(x){var(ABBA_table$fG[pi_cuts == levels(pi_cuts)[x]  & ABBA_table$D >= 0 & ABBA_table$fG >= 0 & ABBA_table$fG <= 1], na.rm = T)})
fhom_vars <- sapply(1:length(levels(pi_cuts)), function(x){var(ABBA_table$fhom[pi_cuts == levels(pi_cuts)[x] & ABBA_table$D >= 0 & ABBA_table$fhom >= 0 & ABBA_table$fhom <= 1], na.rm = T)})

plot(D_vars, type = "b", pch = 19, xlab = expression(paste(pi, " bin")), xaxt = "n", ylab = "Variance", ylim = c(0,0.1))
axis(1, at=1:6, labels = c("0 - 0.01", "0.01 - 0.02", "0.02 - 0.03", "0.03 - 0.04", "0.04 - 0.05", "0.05 - 0.06"), las = 1) 
points(fd_vars, type = "b")
points(fG_vars, type = "b", pch = 2)
points(fhom_vars, type = "b", pch = 3)

mtext(2,text = "A", las = 2, at = 0.1, line = 2.5, cex = 1.5)

legend(3,0.1, legend = c(expression(italic("D")), expression(italic("f"["G"])), expression(italic("f"["hom"])), expression(italic("f"["d"]))), pch = c(19,2,3,1), bty = "n")





ABBA_table <- read.csv("Heliconius_autosome_windows_10kb.csv", as.is = T)
ABBA_table_Z <- read.csv("Heliconius_Zchromosome_windows_10kb.csv", as.is = T)
ABBA_table <- rbind(ABBA_table, ABBA_table_Z)

#filter for windows with enough sites above the minimum data cutoff
ABBA_table <- subset(ABBA_table, sitesOverMinExD >= 6000)


#calculate mean Pi in the three pops

mean_Pi <- apply(ABBA_table[,c("P1_pi", "P2_pi", "P3_pi")], 1, mean, na.rm = T) 

### variances in mean_Pi bins

pi_cuts <- cut(mean_Pi, seq(0,0.04,0.01))

D_vars <- sapply(1:length(levels(pi_cuts)), function(x){var(ABBA_table$D[pi_cuts == levels(pi_cuts)[x] & ABBA_table$D >= 0], na.rm = T)})
fd_vars <- sapply(1:length(levels(pi_cuts)), function(x){var(ABBA_table$fd[pi_cuts == levels(pi_cuts)[x] & ABBA_table$D >= 0] , na.rm = T)})
fG_vars <- sapply(1:length(levels(pi_cuts)), function(x){var(ABBA_table$fG[pi_cuts == levels(pi_cuts)[x]  & ABBA_table$D >= 0 & ABBA_table$fG >= 0 & ABBA_table$fG <= 1], na.rm = T)})
fhom_vars <- sapply(1:length(levels(pi_cuts)), function(x){var(ABBA_table$fhom[pi_cuts == levels(pi_cuts)[x] & ABBA_table$D >= 0 & ABBA_table$fhom >= 0 & ABBA_table$fhom <= 1], na.rm = T)})

plot(D_vars, type = "b", pch = 19, xlab = expression(paste(pi, " bin")), xaxt = "n", ylab = "Variance", ylim = c(0,0.1))
axis(1, at=1:6, labels = c("0 - 0.01", "0.01 - 0.02", "0.02 - 0.03", "0.03 - 0.04", "0.04 - 0.05", "0.05 - 0.06"), las = 1) 
points(fd_vars, type = "b")
points(fG_vars, type = "b", pch = 2)
points(fhom_vars, type = "b", pch = 3)

mtext(2,text = "B", las = 2, at = 0.1, line = 2.5, cex = 1.5)





ABBA_table <- read.csv("Heliconius_autosome_windows_20kb.csv", as.is = T)
ABBA_table_Z <- read.csv("Heliconius_Zchromosome_windows_20kb.csv", as.is = T)
ABBA_table <- rbind(ABBA_table, ABBA_table_Z)

#filter for windows with enough sites above the minimum data cutoff
ABBA_table <- subset(ABBA_table, sitesOverMinExD >= 12000)


#calculate mean Pi in the three pops

mean_Pi <- apply(ABBA_table[,c("P1_pi", "P2_pi", "P3_pi")], 1, mean, na.rm = T) 

### variances in mean_Pi bins

pi_cuts <- cut(mean_Pi, seq(0,0.04,0.01))

D_vars <- sapply(1:length(levels(pi_cuts)), function(x){var(ABBA_table$D[pi_cuts == levels(pi_cuts)[x] & ABBA_table$D >= 0], na.rm = T)})
fd_vars <- sapply(1:length(levels(pi_cuts)), function(x){var(ABBA_table$fd[pi_cuts == levels(pi_cuts)[x] & ABBA_table$D >= 0] , na.rm = T)})
fG_vars <- sapply(1:length(levels(pi_cuts)), function(x){var(ABBA_table$fG[pi_cuts == levels(pi_cuts)[x]  & ABBA_table$D >= 0 & ABBA_table$fG >= 0 & ABBA_table$fG <= 1], na.rm = T)})
fhom_vars <- sapply(1:length(levels(pi_cuts)), function(x){var(ABBA_table$fhom[pi_cuts == levels(pi_cuts)[x] & ABBA_table$D >= 0 & ABBA_table$fhom >= 0 & ABBA_table$fhom <= 1], na.rm = T)})

plot(D_vars, type = "b", pch = 19, xlab = expression(paste(pi, " bin")), xaxt = "n", ylab = "Variance", ylim = c(0,0.1))
axis(1, at=1:6, labels = c("0 - 0.01", "0.01 - 0.02", "0.02 - 0.03", "0.03 - 0.04", "0.04 - 0.05", "0.05 - 0.06"), las = 1) 
points(fd_vars, type = "b")
points(fG_vars, type = "b", pch = 2)
points(fhom_vars, type = "b", pch = 3)

mtext(2,text = "C", las = 2, at = 0.1, line = 2.5, cex = 1.5)






ABBA_table <- read.csv("Heliconius_autosome_windows_50kb.csv", as.is = T)
ABBA_table_Z <- read.csv("Heliconius_Zchromosome_windows_50kb.csv", as.is = T)
ABBA_table <- rbind(ABBA_table, ABBA_table_Z)

#filter for windows with enough sites above the minimum data cutoff
ABBA_table <- subset(ABBA_table, sitesOverMinExD >= 24000)


#calculate mean Pi in the three pops

mean_Pi <- apply(ABBA_table[,c("P1_pi", "P2_pi", "P3_pi")], 1, mean, na.rm = T) 

### variances in mean_Pi bins

pi_cuts <- cut(mean_Pi, seq(0,0.04,0.01))

D_vars <- sapply(1:length(levels(pi_cuts)), function(x){var(ABBA_table$D[pi_cuts == levels(pi_cuts)[x] & ABBA_table$D >= 0], na.rm = T)})
fd_vars <- sapply(1:length(levels(pi_cuts)), function(x){var(ABBA_table$fd[pi_cuts == levels(pi_cuts)[x] & ABBA_table$D >= 0] , na.rm = T)})
fG_vars <- sapply(1:length(levels(pi_cuts)), function(x){var(ABBA_table$fG[pi_cuts == levels(pi_cuts)[x]  & ABBA_table$D >= 0 & ABBA_table$fG >= 0 & ABBA_table$fG <= 1], na.rm = T)})
fhom_vars <- sapply(1:length(levels(pi_cuts)), function(x){var(ABBA_table$fhom[pi_cuts == levels(pi_cuts)[x] & ABBA_table$D >= 0 & ABBA_table$fhom >= 0 & ABBA_table$fhom <= 1], na.rm = T)})

plot(D_vars, type = "b", pch = 19, xlab = expression(paste(pi, " bin")), xaxt = "n", ylab = "Variance", ylim = c(0,0.1))
axis(1, at=1:6, labels = c("0 - 0.01", "0.01 - 0.02", "0.02 - 0.03", "0.03 - 0.04", "0.04 - 0.05", "0.05 - 0.06"), las = 1) 
points(fd_vars, type = "b")
points(fG_vars, type = "b", pch = 2)
points(fhom_vars, type = "b", pch = 3)

mtext(2,text = "D", las = 2, at = 0.1, line = 2.5, cex = 1.5)

dev.off()
