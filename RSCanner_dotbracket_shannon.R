#!/usr/bin/env Rscript
args=commandArgs(TRUE)
if (length(args)<1){
  cat("RSCanner_dotbracket_shannon.R dotbracket.fasta shannon.txt\n")
  cat("\nInput:\n")
  cat("    dotbracket.fasta          - dotbracket file in FASTA format\n")
  cat("    shannon.txt               - single column of shannon entropy values with header 'Shannon'\n")  
  cat("\nOutput:\n")
  cat("    heatmap.tiff              - output heatmap figure\n")
  quit()
}

library(readxl)
library(MASS)
library(tidyverse)
library(dplyr)
library(grid)
library(boxplotdbl)
library(xlsx)
library(seqinr)

#USER: Input dot-bracket file
dotbracketinput <- seqinr::read.fasta(args[1], set.attributes = FALSE, whole.header = FALSE)
dotbracket_vector <- dotbracketinput$ENERGY[((length(dotbracketinput$ENERGY)/2)+1):length(dotbracketinput$ENERGY)]
#USER: Input ct file
#ctfileinput <- read.table("HCV_JFH1_2a_genome_model.ct")
#USER (OPTIONAL): Input list of Shannon values
#shannoninput <- (read_excel("Dataset_S1_HCV_genome_SHAPE_reactivities.xlsx", sheet = 3, col_names = TRUE))
shannoninput <- read_table(args[1], header = TRUE)
shannon <- shannoninput$Shannon
#USER: Input window size for BPC calculation
window_size <- as.integer(readline(prompt = "Input window size for BPC calculation: "))
#USER: Input window size for Shannon entropy smoothing
window_size_shan <- as.integer(readline(prompt = "Input window size for Shannon Entropy smoothing: "))
#USER: Input Shannon Entropy cutoff
SE_cutoff <- as.double(readline(prompt = "Input Shannon Entropy cutoff: "))
#USER: Input BPC cutoff
BPC_cutoff <- as.double(readline(prompt = "Input Shannon Entropy cutoff: "))

dotbracket_vector_bin <- numeric(length(dotbracket_vector))
for (i in 1:length(dotbracket_vector)) {
  if (dotbracket_vector[i] == ".") {
    dotbracket_vector_bin[i] <- 1
  } else {
    dotbracket_vector_bin[i] <- 0
  }
}
if ((window_size %% 2) == 0) {
  gapnum <- window_size/2 + 1
} else {
  gapnum <- (window_size+1)/2 + 1
}
window_initvec <- numeric(length(dotbracket_vector))
window_termvec <- numeric(length(dotbracket_vector))
window_termvec[1] <- gapnum
for (j in 1:gapnum) {
  window_initvec[j] <- 1
}
for (k in (gapnum + 1):(length(dotbracket_vector))){
  window_initvec[k] <- window_initvec[k-1] + 1
}
for (l in (length(dotbracket_vector) - gapnum + 1):(length(dotbracket_vector))) {
  window_termvec[l] <- length(dotbracket_vector)
}
for (m in 2:(length(dotbracket_vector) - gapnum)) {
  window_termvec[m] <- window_termvec[m-1] + 1
}
dotcount <- numeric(length(dotbracket_vector))
for (o in 1:length(dotbracket_vector)) {
  dotcount[o] <- sum(dotbracket_vector_bin[window_initvec[o]:window_termvec[o]])
}
dotperc <- numeric(length(dotbracket_vector))
for (p in 1:length(dotbracket_vector)) {
  dotperc[p] <- (dotcount[p])/(window_termvec[p] - window_initvec[p] + 1)
}
oneminus_dotperc <- 1 - dotperc

#Shannon smoothing with customized size sliding window
#movmedian calculation
if ((window_size_shan %% 2) == 0) {
  gapnum_shan <- window_size_shan/2 + 1
} else {
  gapnum_shan <- (window_size_shan+1)/2 + 1
}
window_initvec_shan <- numeric(length(shannon))
window_termvec_shan <- numeric(length(shannon))
window_termvec_shan[1] <- gapnum_shan
for (j in 1:gapnum_shan) {
  window_initvec_shan[j] <- 1
}
for (k in (gapnum_shan + 1):(length(shannon))){
  window_initvec_shan[k] <- window_initvec_shan[k-1] + 1
}
for (l in (length(shannon) - gapnum_shan + 1):(length(shannon))) {
  window_termvec_shan[l] <- length(shannon)
}
for (m in 2:(length(shannon) - gapnum_shan)) {
  window_termvec_shan[m] <- window_termvec_shan[m-1] + 1
}
med_shan <- numeric(length(shannon))
for (o in 1:length(shannon)) {
  med_shan[o] <- median(shannon[window_initvec_shan[o]:window_termvec_shan[o]])
}

a <- which(oneminus_dotperc > quantile(oneminus_dotperc, probs=c(BPC_cutoff),name=FALSE)) #indices of the bpcs that are above cutoff
b <- which(med_shan < quantile(med_shan, probs=c(SE_cutoff),name=FALSE)) #indices of shannons that are below cutoff
abinter <- intersect(a,b) #intersection of the two vectors

#plot BPC along the full length RNA 
bpcplot <- plot(oneminus_dotperc, type = 'l')
ggsave(bpcplot, device="tiff", width=7, height=3, dpi=300)

#plot the smoothed Shannon entropy
shannonplot <- plot(med_shan, type = 'l')
ggsave(shannonplot, device="tiff", width=7, height=3, dpi=300)

#plot the percentage of well-defined structures in non-overlaping bins along the RNA
finalwind <- 100
finalwind_init <- numeric(ceiling(length(shannon)/finalwind))
finalwind_fin <- numeric(ceiling(length(shannon)/finalwind))
finalwind_init[1] <- 1
for (k in 2:length(finalwind_init)) {
  finalwind_init[k] <- finalwind_init[k-1] + finalwind
}
finalwind_fin[1] <- finalwind
for (k in 2:length(finalwind_fin)) {
  finalwind_fin[k] <- finalwind_fin[k-1] + finalwind
}
finalwind_fin[length(finalwind_fin)] <- length(shannon)
structure_counts <- numeric(length(finalwind_fin))
for (o in 1:length(structure_counts)) {
  count_raw <- sum(finalwind_init[o]<=abinter & abinter<=finalwind_fin[o])
  window_size <- (finalwind_fin[o] - finalwind_init[o] + 1)
  structure_counts[o] <- 100*(count_raw/window_size)
}
finalwind_inds <- seq(from = finalwind/2, to = length(shannon), by = finalwind)

if (length(finalwind_inds) != length(structure_counts)) {
  finalwind_inds_real <- numeric(length(finalwind_inds)+1)
  finalwind_inds_real[1:length(finalwind_inds_real)-1] <- finalwind_inds
  finalwind_inds_real[length(finalwind_inds_real)] <- finalwind_inds[length(finalwind_inds)]+finalwind
} else (finalwind_inds_real <- finalwind_inds)

#color ramp creation
colorramp <-  colorRampPalette(colors=c("#FFFF00", "#FF0000"))(101)
inds_colors <- numeric(length(finalwind_inds))
for (i in 1:length(structure_counts)) {
  inds_colors[i] <- colorramp[structure_counts[i]+1]
}
inds_colors[which(inds_colors == "#FFFF00")] <- "#FFFFFF"
dat <- data.frame(pos = finalwind_inds, vals = structure_counts, cols = inds_colors)
start <- finalwind_inds - 50
end <- finalwind_inds + 50
end[length(end)] <- length(shannon)
## highlight region data
rects <- data.frame(start=start, end=end, group=seq_along(start))
# min and max x and y values
ymin <- min(dat$vals)
ymax <- max(dat$vals)
xmin <- start[1]
xmax <- end[length(end)]
library(ggplot2)
heatmap <- ggplot(data=dat, aes(pos, vals)) +
  theme_classic() +
  geom_rect(data=rects, inherit.aes=FALSE, aes(xmin=start, xmax=end, ymin=ymin,
                                               ymax=ymax, group=group), color="transparent", 
            fill=as.character(dat$cols), alpha=1) + geom_line(lty=1, color="black", size = 0.9) +
  coord_cartesian(xlim = c(xmin, xmax+50), ylim = c(ymin, ymax)) +
  xlab("Nucleotide") + ylab("% Structure Content") + theme(axis.text = element_text(size = 10, color="black"), 
                                                           axis.title = element_text(size = 12), panel.border = element_rect(color="black", fill=NA, size = 1))

ggsave(heatmap, device="tiff", width=7, height=3, dpi=300)