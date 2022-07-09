#!/usr/bin/env Rscript
args=commandArgs(TRUE)
if (length(args)<1){
  cat("To run the script, do\n")
  cat("\n")
  cat("    Rscript RSCanner_CT_shannon.R <input_structure.ct>  <shannon.txt>\n")
  cat("\nInput:\n")
  cat("    input_structure.ct        - CT format file; six columns after one header line\n")
  cat("    shannon.txt               - two columns in a tab delimited file, col1 = index, col2 = shannon entropy values with no header \n")  
  cat("\nOutput:\n")
  cat("    bpcplot.tiff                             - output BPC line plot figure\n")
  cat("    smoothed_Shannonplot.tiff                - output smoothed SE line plot figure\n")
  cat("    structure_counts_histogram.tiff          - output histogram figure\n")
  cat("    structure_counts_heatmap.tiff            - output final line plot/heatmap figure\n")
  cat("\n")
  cat("    bpc_data.csv                             - output tab delimited file, col1=index, col2=nucleotide number, col3=BPC\n")
  cat("    smoothed_Shannon_data.csv                - output tab delimited file, col1=index, col2=nucleotide number, col3=smoothed Shannon Entropy\n")
  cat("    structure_counts.csv                     - output tab delimited file, col1=index, col2=structure counts, col3=BPC, col4=smoothed Shannon Entropy\n")
  cat("    ordered_structure_table.csv              - output tab delimited file, col1=index, col2=bin number, col3=% structure content, col4=bin start (nt), col5=bin end (nt)\n")
  quit()
}

#Load the required packages
library(tidyverse)
library(dplyr)
library(ggplot2)

## Step 1: Inputs

#USER: Input ct file
ctfileinput <- read.table(args[1], skip = 1)$V5

#USER (OPTIONAL): Input list of Shannon values
shannoninput <- read.table(args[2], header = FALSE)
shannon <- shannoninput$V2

#USER: Input window size for BPC calculation
cat("Input integer window size for BPC calculation (use 51 as default): ")
window_size <- as.integer(readLines("stdin", 1))
#USER: Input window size for Shannon entropy smoothing
cat("Input integer window size for Shannon Entropy smoothing (use 51 as default): ")
window_size_shan <- as.integer(readLines("stdin", 1))
#USER: Input Shannon Entropy cutoff
cat("Input Shannon Entropy percentile cutoff, decimal between 0 and 1 (use 0.5 as default): ")
SE_cutoff <- as.double(readLines("stdin", 1))
#USER: Input BPC cutoff
cat("Input BPC percentile cutoff, decimal between 0 and 1 (use 0.5 as default): ")
BPC_cutoff <- as.double(readLines("stdin", 1))

#Convert input ct file to dotbracket vector format
dotbracket_nums <- ctfileinput
indexnotzero <- which(dotbracket_nums != 0)
for (i in 2:length(indexnotzero)) {
  if(indexnotzero[i] - indexnotzero[i-1] != 1) {
    dotbracket_nums[indexnotzero[i]] <- "*"
  }
  dotbracket_nums[indexnotzero[1]] <- "*"
}
wichaster <- which(dotbracket_nums == "*")
for (j in 1:length(wichaster)) {
  if (j %% 2 == 0) {
    dotbracket_nums[wichaster[j]] <- ")"
  } else {
    dotbracket_nums[wichaster[j]] <- "("
  }
}
dotbracket_real <- dotbracket_nums
`%!in%` <- Negate(`%in%`)
for (k in 1:length(dotbracket_nums)) {
  ind <- wichaster[which(wichaster < k)[length(which(wichaster < k))]]
  if (dotbracket_real[k] %!in% c("0", "(", ")")) {
    dotbracket_real[k] <- dotbracket_nums[ind]
  } else {
    dotbracket_real[k] <- dotbracket_nums[k]
  }
}
dotbracket_real[which(dotbracket_real == "0")] <- "."

dotbracket_vector <- dotbracket_real

## Step 2: Sliding window calculations

#Base pair content (BPC) calculation with customized sliding window size
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
  gapnum <- (window_size+1)/2 
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

#Shannon entropy (SE) smoothing with customized sliding window size
if ((window_size_shan %% 2) == 0) {
  gapnum_shan <- window_size_shan/2 + 1
} else {
  gapnum_shan <- (window_size_shan+1)/2 
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

## Step 3: Filtering

#Calculate proportion of regions that have BPC above the threshold and smoothed SE below the threshold ("structure counts")
a <- which(oneminus_dotperc > quantile(oneminus_dotperc, probs=c(BPC_cutoff),name=FALSE)) #indices of the bpcs that are above cutoff
b <- which(med_shan < quantile(med_shan, probs=c(SE_cutoff),name=FALSE)) #indices of shannons that are below cutoff
abinter <- intersect(a,b) #intersection of the two vectors

passing_nucs <- as.data.frame(cbind(`Nucleotide` = abinter, 
                                    `Base Pair Content` = oneminus_dotperc[abinter], 
                                    `Shannon` = med_shan[abinter]))

#Save csv of nucleotides that pass both thresholds
write.csv(passing_nucs, "structure_counts.csv")

#USER: input the x-axis bounds (nt) for visualization
cat("The transcript inputted is of length", length(shannon), "\n\n")
cat("Please input desired bounds for visualization, between 1 and", length(shannon), "\n\n")
cat("Input integer lower bound (default value of 1): ")
x_low_bound <- as.integer(readLines("stdin", 1))
cat("Input integer upper bound (default value of", length(shannon), "): ")
x_upper_bound <- as.integer(readLines("stdin", 1))

if (x_upper_bound > length(shannon)) {
  x_upper_bound <- length(shannon)
} else {x_upper_bound <- x_upper_bound}
pdf(NULL)

#Plot BPC along the full length of the RNA 
bpcplotdata <- as.data.frame(cbind(seq(1, length(oneminus_dotperc)), oneminus_dotperc)) %>% rename("Nucleotide Position" = V1) %>% 
  rename("Base Pair Content" = oneminus_dotperc)
ggplot(data = bpcplotdata, aes(x = `Nucleotide Position`, y = `Base Pair Content`)) + theme_classic() + 
  geom_line() + 
  geom_hline(yintercept = quantile(oneminus_dotperc, probs=c(BPC_cutoff),name=FALSE), linetype = "dashed", color = "black", size = 1) + 
  scale_x_continuous(limits = c(x_low_bound, x_upper_bound))+
  theme(axis.text = element_text(size = 10, color="black"), axis.title = element_text(size = 12), panel.border = element_rect(color="black", fill=NA, size = 1))
pdf(NULL)

#Save BPC data as a csv
write.csv(bpcplotdata, "bpc_data.csv")

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

#Save BPC lineplot
ggsave("bpcplot.tiff", device="tiff", width=widthinput, height=heightinput, dpi=dpiinput)

#Plot smoothed SE along the full length of the RNA
shanplotdata <- as.data.frame(cbind(seq(1, length(med_shan)), med_shan)) %>% rename("Nucleotide Position" = V1) %>% 
  rename("Smoothed Median Shannon Entropy" = med_shan)
ggplot(data = shanplotdata, aes(x = `Nucleotide Position`, y = `Smoothed Median Shannon Entropy`)) + theme_classic() +
  geom_line() + 
  geom_hline(yintercept = quantile(med_shan, probs=c(SE_cutoff),name=FALSE), linetype = "dashed", color = "black", size = 1) + 
  scale_x_continuous(limits = c(x_low_bound, x_upper_bound))+
  theme(axis.text = element_text(size = 10, color="black"), axis.title = element_text(size = 12), panel.border = element_rect(color="black", fill=NA, size = 1))

#Save smoothed SE data as a csv
write.csv(shanplotdata, "smoothed_Shannon_data.csv")

cat("\n Computation complete... saving shannon image. \n\n")

#Save smoothed SE lineplot
ggsave("smoothed_Shannonplot.tiff", device="tiff", width=widthinput, height=heightinput, dpi=dpiinput)

## Step 4: Binning and Visualization

#USER: Input window size for final computation
cat("Input integer window length for histogram and heatmap (use 100 as default): ")
finalwind <- as.integer(readLines("stdin", 1))

if ((finalwind %% 2) == 1) {
  finalwind <- finalwind + 1
} else {
  finalwind <- finalwind
}

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

#Save the ordered table to a csv
write.csv(ordered_results_table, "ordered_structure_table.csv")

#Plot the percentage of structure content in non-overlapping bins along the RNA and overlay with a heatmap
half_finalwind <- finalwind/2

#Color ramp creation
colorramp <-  colorRampPalette(colors=c("#FFFF00", "#FF0000"))(finalwind+1)
inds_colors <- numeric(length(finalwind_inds_real))
for (i in 1:length(structure_counts)) {
  inds_colors[i] <- colorramp[structure_counts[i]+1]
}
inds_colors[which(inds_colors == "#FFFF00")] <- "#FFFFFF"
dat <- data.frame(pos = finalwind_inds_real, vals = structure_counts, cols = inds_colors)
starty <- finalwind_inds_real - half_finalwind
endy <- finalwind_inds_real + half_finalwind
endy[length(endy)] <- length(shannon)
#Highlight region data
rects <- data.frame(start=starty, end=endy, group=seq_along(starty))

rects_lowwy <- which(rects$start == floor(x_low_bound/finalwind)*100)

if ((ceiling(x_upper_bound/finalwind)*100) > rects[nrow(rects),]$end) {
  rects_uppy <- which(rects$end == rects[nrow(rects),]$end)
} else {
  rects_uppy <- which(rects$end == ceiling(x_upper_bound/finalwind)*100)
}

new_rects <- rects %>% slice(rects_lowwy:rects_uppy)

new_dat <- dat %>% slice(rects_lowwy:rects_uppy)

#Min and max x and y values
ymin <- min(dat$vals)
ymax <- max(dat$vals)
xmin <- new_rects$start[1]
xmax <- new_rects$end[length(new_rects$end)]

#Generate line plot with heatmap overlay
heatmap <- ggplot(data=new_dat, aes(pos, vals)) +
  theme_classic() +
  geom_rect(data=new_rects, inherit.aes=FALSE, aes(xmin=starty[rects_lowwy:rects_uppy], xmax=endy[rects_lowwy:rects_uppy], ymin=ymin,
                                                   ymax=ymax, group=group), color="transparent", 
            fill=as.character(new_dat$cols), alpha=1) + geom_line(lty=1, color="black", size = 0.9) +
  coord_cartesian(xlim = c(xmin, xmax+half_finalwind), ylim = c(ymin, ymax)) +
  xlab("Nucleotide") + ylab("% Structure Content") + theme(axis.text = element_text(size = 10, color="black"), 
                                                           axis.title = element_text(size = 12), panel.border = element_rect(color="black", fill=NA, size = 1))

cat("\n Computation complete... saving heatmap image. \n\n")

#Save line plot/heatmap
ggsave("structure_counts_heatmap.tiff", device="tiff", width=widthinput, height=heightinput, dpi=dpiinput)

#Plot histogram of structure counts along the RNA
bar <- ggplot(data=(new_dat %>% rename(`Structure Counts` = vals) %>% rename(`Nucleotide` = pos)), aes(x=`Nucleotide`, y=`Structure Counts`)) +
  geom_bar(stat="identity",  fill="grey", colour="black", width=100)+
  theme_classic()+
  theme(axis.text = element_text(size = 10, color="black"), axis.title = element_text(size = 12), panel.border = element_rect(color="black", fill=NA, size = 1))

#Save histogram
ggsave("structure_counts_histogram.tiff", device="tiff", width=widthinput, height=heightinput, dpi=dpiinput)
