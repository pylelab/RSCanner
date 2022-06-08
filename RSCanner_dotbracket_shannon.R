#!/usr/bin/env Rscript
args=commandArgs(TRUE)
if (length(args)<1){
  cat("RSCanner_dotbracket_shannon.R dotbracket.fasta shannon.txt\n")
  cat("\nInput:\n")
  cat("    dotbracket.fasta          - dotbracket file in FASTA format; three lines: title, sequence, dotbracket structure\n")
  cat("    shannon.txt               - two columns in a tab delimited file, col1 = index, col2 = shannon entropy values with no header \n")  
  cat("\nOutput:\n")
  cat("    bpcplot.tiff                             - output heatmap figure\n")
  cat("    shannonplot.tiff                         - output heatmap figure\n")
  cat("    heatmap.tiff                             - output heatmap figure\n")
  cat("    ordered_structure_table.csv              - output heatmap figure\n")
  quit()
}

library(tidyverse)
library(dplyr)

#USER: Input dot-bracket file
#dotbracket file must have three lines: title line, sequence line, dotbracket line, and no newline characters except 
#for the ones separating these three lines.
dotbracketinput <- readLines(args[1])[3]
dotbracket_vector <- strsplit(dotbracketinput, split = NULL)[[1]]

#USER (OPTIONAL): Input list of Shannon values
shannoninput <- read.table(args[2], header = FALSE)
shannon <- shannoninput$V2

#USER: Input window size for BPC calculation
cat("Input integer window size for BPC calculation (use 50 as default): ")
window_size <- as.integer(readLines("stdin", 1))
#USER: Input window size for Shannon entropy smoothing
cat("Input integer window size for Shannon Entropy smoothing (use 50 as default): ")
window_size_shan <- as.integer(readLines("stdin", 1))
#USER: Input Shannon Entropy cutoff
cat("Input Shannon Entropy percentile cutoff, decimal between 0 and 1 (use 0.5 as default): ")
SE_cutoff <- as.double(readLines("stdin", 1))
#USER: Input BPC cutoff
cat("Input BPC percentile cutoff, decimal between 0 and 1 (use 0.5 as default): ")
BPC_cutoff <- as.double(readLines("stdin", 1))

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
bpcplotdata <- as.data.frame(cbind(seq(1, length(oneminus_dotperc)), oneminus_dotperc)) %>% rename("Nucleotide Position" = V1) %>% 
  rename("Base Pair Content" = oneminus_dotperc)
bpcplot <- ggplot(data = bpcplotdata, aes(x = `Nucleotide Position`, y = `Base Pair Content`)) + theme_classic() +
  geom_line() + geom_hline(yintercept = quantile(oneminus_dotperc, probs=c(BPC_cutoff),name=FALSE), linetype = "dashed", color = "blue", size = 1)

cat("\n Input image save settings: \n\n")
#USER: Input width of saved image
cat("Input integer width (use 7 as default): ")
widthinput <- as.integer(readLines("stdin", 1))
#USER: Input height of saved image
cat("Input integer height (use 3 as default): ")
heightinput <- as.integer(readLines("stdin", 1))
#USER: Input resolution of saved image
cat("Input integer resolution (dpi) (use 300 as default): ")
dpiinput <- as.integer(readLines("stdin", 1))

cat("\n Computation complete... saving base-pair content image. \n\n")

ggsave("bpcplot.tiff", device="tiff", width=widthinput, height=heightinput, dpi=dpiinput)

#plot the smoothed Shannon entropy
shanplotdata <- as.data.frame(cbind(seq(1, length(med_shan)), med_shan)) %>% rename("Nucleotide Position" = V1) %>% 
  rename("Smoothed Median Shannon Entropy" = med_shan)
shannonplot <- ggplot(data = shanplotdata, aes(x = `Nucleotide Position`, y = `Smoothed Median Shannon Entropy`)) + theme_classic() +
  geom_line() + geom_hline(yintercept = quantile(med_shan, probs=c(SE_cutoff),name=FALSE), linetype = "dashed", color = "blue", size = 1)

cat("\n Computation complete... saving shannon image. \n\n")

ggsave("shannonplot.tiff", device="tiff", width=widthinput, height=heightinput, dpi=dpiinput)

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
} else {
  finalwind_inds_real <- finalwind_inds
  }

bin_number <- seq(from = 1, to = length(finalwind_inds_real), by = 1)
unordered_results_table <- as.data.frame(cbind(bin_number, structure_counts, finalwind_init, finalwind_fin))
ordered_results_table <- arrange(unordered_results_table, desc(structure_counts)) %>% rename("Bin Number" = bin_number) %>%
  rename("% Structure Content" = structure_counts) %>% rename("Bin Start (nt)" = finalwind_init) %>% 
  rename("Bin End (nt)" = finalwind_fin)

write.csv(ordered_results_table, "ordered_structure_table.csv")

#color ramp creation
colorramp <-  colorRampPalette(colors=c("#FFFF00", "#FF0000"))(101)
inds_colors <- numeric(length(finalwind_inds_real))
for (i in 1:length(structure_counts)) {
  inds_colors[i] <- colorramp[structure_counts[i]+1]
}
inds_colors[which(inds_colors == "#FFFF00")] <- "#FFFFFF"
dat <- data.frame(pos = finalwind_inds_real, vals = structure_counts, cols = inds_colors)
start <- finalwind_inds_real - 50
end <- finalwind_inds_real + 50
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

cat("\n Computation complete... saving heatmap image. \n\n")

ggsave("heatmap.tiff", device="tiff", width=widthinput, height=heightinput, dpi=dpiinput)
