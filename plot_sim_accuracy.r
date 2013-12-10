#!/usr/bin/env Rscript
# plot_sim_accuracy.r
# Plot accuracy for all models with oldest split at 2MYA

# Written for "Evaluating statistics for the identification of introgressed loci"
# by Simon H. Martin, John W. Davey and Chris D. Jiggins
# John Davey:   jd626@cam.ac.uk
# Simon Martin: shm45@cam.ac.uk
# November-December 2013


models<-read.delim("model_files_win10000_size5000_ne1000000_2013-11-29.partition.summary.tsv")


mf_models <- subset(models, Stat == "mfD0")

D_models <- subset(models, Stat == "D")



### plots using two most informative variables

#we fix one parameter, the deepest split. We plot all available models with Alt23 on the X axis, and Bgt12 a different lines

#vector of gradient colour names
gradCols <- vector(length = 10)
splits <- round(seq(0.2,2,0.2),1)
for (i in 1:length(splits))  gradCols[i] <- rgb(1-(splits[i]/2),1-(splits[i]/2),1-(splits[i]/2),1)


pdf(file = "model_accuracy_plots.pdf", width = 7, height = 10)

par(mfrow = c(3,2), mar = c(5.5,4,1.5,1))

#gf23

D_models_Bt123_2_GF23 <- subset(D_models, Background_t123 == 2 & ModelType == "gene flow P2-P3")
plot(0,0,cex = 0, xlim = c(0,2), ylim = c(0,100), xlab = "Alternate t23", ylab = "percent predicted by D", main = "gene flow P2->P3 (Background t123 = 2)", bty = "n")
for (i in c(0.4,0.6,0.8,1,1.2,1.4,1.6,1.8)){ lines(D_models_Bt123_2_GF23$Alternate_t23[D_models_Bt123_2_GF23$Background_t21 == i], D_models_Bt123_2_GF23$OutlierAlternatePC[D_models_Bt123_2_GF23$Background_t21 == i], type = "b", pch = 21, col = gradCols[splits==i], bg = gradCols[splits==i], cex = 1.8) }
legend(1.4,100,legend = splits[2:9], pch = 19, cex = 1.3, col = gradCols[2:9], title = "Background t12", bty = "n", ncol = 2)

mf_models_Bt123_2_GF23 <- subset(mf_models, Background_t123 == 2 & ModelType == "gene flow P2-P3")
plot(0,0,cex = 0, xlim = c(0,2), ylim = c(0,100), xlab = "Alternate t23", ylab = "percent predicted by mf", main = "gene flow P2->P3 (Background t123 = 2)", bty = "n")
for (i in c(0.4,0.6,0.8,1,1.2,1.4,1.6,1.8)){  lines(mf_models_Bt123_2_GF23$Alternate_t23[mf_models_Bt123_2_GF23$Background_t21 == i], mf_models_Bt123_2_GF23$OutlierAlternatePC[mf_models_Bt123_2_GF23$Background_t21 == i], type = "b", pch = 21, col = gradCols[splits==i], bg = gradCols[splits==i], cex = 1.8)  }
legend(1.4,100,legend = splits[2:9], pch = 19, cex = 1.3, col = gradCols[2:9], title = "Background t12", bty = "n", ncol = 2)


#gf32

D_models_Bt123_2_GF32 <- subset(D_models, Background_t123 == 2 & ModelType == "gene flow P3-P2")
plot(0,0,cex = 0, xlim = c(0,2), ylim = c(0,100), xlab = "Alternate t23", ylab = "percent predicted by D", main = "gene flow P3->P2 (Background t123 = 2)", bty = "n")
for (i in c(0.4,0.6,0.8,1,1.2,1.4,1.6,1.8)){ lines(D_models_Bt123_2_GF32$Alternate_t23[D_models_Bt123_2_GF32$Background_t21 == i], D_models_Bt123_2_GF32$OutlierAlternatePC[D_models_Bt123_2_GF32$Background_t21 == i], type = "b", pch = 21, col = gradCols[splits==i], bg = gradCols[splits==i], cex = 1.8) }
legend(1.4,100,legend = splits[2:9], pch = 19, cex = 1.3, col = gradCols[2:9], title = "Background t12", bty = "n", ncol = 2)

mf_models_Bt123_2_GF32 <- subset(mf_models, Background_t123 == 2 & ModelType == "gene flow P3-P2")
plot(0,0,cex = 0, xlim = c(0,2), ylim = c(0,100), xlab = "Alternate t23", ylab = "percent predicted by mf", main = "gene flow P3->P2 (Background t123 = 2)", bty = "n")
for (i in c(0.4,0.6,0.8,1,1.2,1.4,1.6,1.8)){  lines(mf_models_Bt123_2_GF32$Alternate_t23[mf_models_Bt123_2_GF32$Background_t21 == i], mf_models_Bt123_2_GF32$OutlierAlternatePC[mf_models_Bt123_2_GF32$Background_t21 == i], type = "b", pch = 21, col = gradCols[splits==i], bg = gradCols[splits==i], cex = 1.8)  }
legend(1.4,100,legend = splits[2:9], pch = 19, cex = 1.3, col = gradCols[2:9], title = "Background t12", bty = "n", ncol = 2)


#structure

D_models_At123_2_AS <- subset(D_models, Alternate_t123 == 2 & ModelType == "Ancestral structure")
plot(0,0,cex = 0, xlim = c(0,2), ylim = c(0,100), xlab = "Alternate t23", ylab = "percent predicted by D", main = "ancestral structure (Alternate t123 = 2)", bty = "n")
for (i in c(0.2,0.4,0.6,0.8,1,1.2,1.4,1.6)){ lines(D_models_At123_2_AS$Alternate_t23[D_models_At123_2_AS$Background_t21 == i], D_models_At123_2_AS$OutlierAlternatePC[D_models_At123_2_AS$Background_t21 == i], type = "b", pch = 21, col = gradCols[which(splits==i)+1], bg = gradCols[which(splits==i)+1], cex = 1.8) }
legend(1.4,100,legend = splits[1:8], pch = 19, cex = 1.3, col = gradCols[2:9], title = "Background t12", bty = "n", ncol = 2)

mf_models_At123_2_AS <- subset(mf_models, Alternate_t123 == 2 & ModelType == "Ancestral structure")
plot(0,0,cex = 0, xlim = c(0,2), ylim = c(0,100), xlab = "Alternate t23", ylab = "percent predicted by mf", main = "ancestral structure (Alternate t123 = 2)", bty = "n")
for (i in c(0.2,0.4,0.6,0.8,1,1.2,1.4,1.6)){  lines(mf_models_At123_2_AS$Alternate_t23[mf_models_At123_2_AS$Background_t21 == i], mf_models_At123_2_AS$OutlierAlternatePC[mf_models_At123_2_AS$Background_t21 == i], type = "b", pch = 21, col = gradCols[which(splits==i)+1], bg = gradCols[which(splits==i)+1], cex = 1.8)  }
legend(1.4,100,legend = splits[1:8], pch = 19, cex = 1.3, col = gradCols[2:9], title = "Background t12", bty = "n", ncol = 2)


dev.off()

