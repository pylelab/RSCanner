# RSCanner

RSCanner is an R package for rapid assessment and visualization of RNA structural content that is particularly useful for long RNAs. The script has been written as an Rmd file, and command line functionality has been provided via two separate R scripts, to be used depending on the type of input given. You can access fully functional RSCanner using either 1) the R scripts via a Linux terminal (recommended for new R users) or 2) the Rmd notebook via RStudio (a popular R IDE). See usage instructions below. 

## Inputs and algorithm description

Input files: 1 - RNA secondary structure in either dot bracket notation (FASTA file) or CT (connectivity table) format; 2 - positional Shannon entropy for each nucleotide in the RNA, calculated from base pairing probabilities.

RSCanner scans the secondary structure provided along with Shannon entropy values in sliding windows covering the entire RNA and calculates the base pair content (BPC) and median entropy values for each window. Then, it computes nucleotide positions with BPC values above a user-defined cutoff (default = 50th percentile) and Shannon entropy values below a user-defined cutoff (default = 50th percentile). These positions are termed structure counts. RSCanner then outputs the frequency distribution of structure counts across the RNA as both heatmap and histogram plots.  

RSCanner uses a centered sliding window (default=51nt, 25nt flanks to each side of the center nucleotide, moved in steps of 1 nt across the RNA) to calculate BPC and median Shannon Entropy. If an even number window size is inputted, RSCanner will use the next odd number as the window size (i.e., if input=50, the actual window size=51). In this way, a centered sliding window is always used in all calculations. When dealing with the ends, RSCanner will gradually truncate the window size (i.e., for 5' end windows - window #26=[1nt,51nt], window#25=[1nt,50nt], window#24 =[1nt,49nt], etc; 3’ end windows are computed analogously).

## Installation
RSCanner requires pre-installation of the R computing language and three R packages.

### Step 1: Install R computing language
From the R Project documentation, "R is a free software environment for statistical computing and graphics. It compiles and runs on a wide variety of UNIX platforms, Windows and MacOS." 

To use our RSCanner tool, download R from https://www.r-project.org/ by choosing a CRAN mirror. For example, you can use the CRAN mirror hosted by the National Institute for Computational Sciences at https://mirrors.nics.utk.edu/cran/. Choose your machine and operating system, and download R following the directions on the documentation.

### Step 2: Clone the RSCanner Github repository by the Pyle Lab at Yale University
Clone this entire Github repository to your local machine. The repository will download as a zip file called 'RSCanner-main.zip' by default. Unzip the file. We will proceed with usage instructions assuming the user is in 'RSCanner-main' directory on their machine.

### (Not required) Step 3: Install R packages

#### The Rscripts automatically check for the required packages and install them if needed. However, we included instructions below for manual installation.

The installation instructions for the R packages differ depending on whether you are using the R scripts or the Rmd notebook in RStudio. The RSCanner script has three R package dependencies: 'tidyverse', 'dplyr', and 'ggplot2'. The user can install these packages by opening the R Console application that they installed in Step 1 above, and running the following lines on the R console command line that opens.

#### Step 3.1: Open R
After clicking the R application that you installed in Step 1 above, you will see an *R Console* as follows:
```
>
```
#### Step 3.2: Install Package
For example, run:
```
> install.packages("tidyverse")
```

#### Step 3.3: Select CRAN Mirror
The R Console will open a selection menu, from which the user must select a CRAN mirror (any will do, the most convenient is to select the CRAN mirror corresponding to your geographical area, e.g. USA)

#### Step 3.4: Repeat steps 3.1-3.3 for packages 'dplyr' and 'ggplot2' as well.
i.e. run
```
> install.packages("dplyr")
> install.packages("ggplot2")
```
and then select mirror.

We're now ready to go! The R language and all package dependencies have been installed. Now, you can use RSCanner from your terminal.
Navigate to your RSCanner-main directory that you have cloned in Step 2 above and proceed as follows.

## Usage: RSCanner_CT_shannon.R
To be used when the secondary structure input is a CT file.
For usage information and detailed input/output information, run the following in your terminal:

```
 Rscript RSCanner_CT_shannon.R
 
 # Terminal output
 
 To run the script, do

     Rscript RSCanner_CT_shannon.R <input_structure.ct>  <shannon.txt>

 Input:
     input_structure.ct        - CT format file; six columns after one header line
     shannon.txt               - two columns in a tab delimited file, col1 = index, col2 = shannon entropy values with no header 

 Output:
     bpcplot.tiff                             - output BPC line plot figure, dashed line represents the global median
     smoothed_Shannonplot.tiff                - output smoothed SE line plot figure, dashed line represents the global median
     structure_counts_histogram.tiff          - output histogram figure
     structure_counts_heatmap.tiff            - output final line plot/heatmap figure

     bpc_data.csv                             - output tab delimited file, col1=index, col2=nucleotide number, col3=BPC
     smoothed_Shannon_data.csv                - output tab delimited file, col1=index, col2=nucleotide number, col3=smoothed Shannon Entropy
     structure_counts.csv                     - output tab delimited file, col1=index, col2=structure counts, col3=BPC, col4=smoothed Shannon Entropy
     ordered_structure_table.csv              - output tab delimited file, col1=index, col2=bin number, col3=% structure content, col4=bin start (nt), col5=bin end (nt)

```

## Usage: RSCanner_dotbracket_shannon.R
To be used when the secondary structure input is in dot bracket format (FASTA file).
For usage information and detailed input/output information, run the following in your terminal:

```
 Rscript RSCanner_dotbracket_shannon.R
```

## User-defined parameters
The program will prompt the user for specific parameters. The recommended values are displayed as the program runs and are listed below:

Integer window size for BPC calculation: 51<br/>
Integer window size for Shannon Entropy smoothing: 51<br/>
Shannon Entropy percentile cutoff, decimal between 0 and 1: 0.5<br/>
BPC percentile cutoff, decimal between 0 and 1: 0.5<br/>
Integer lower bound: 1<br/>
Integer upper bound: length of transcript<br/>
Integer width of images: 7<br/>
Integer height of images: 3<br/>
Integer resolution of images (dpi): 300<br/>
Integer window length for histogram and heatmap: 100<br/>


## Example
Here is a fully worked example using the sample CT, dotbracket, and shannon text files included in this repository. Simply download the entire repository as-is, and run the following commands from within this repo directory.

### Using RSCanner_CT_shannon.R 
```
mkdir sample_outputs

cd sample_outputs

Rscript ../RSCanner_CT_shannon.R ../sample_test_files/HCV_Jc1.ct ../sample_test_files/HCV_shannon.txt

```

### Using RSCanner_dotbracket_shannon.R
```
mkdir sample_outputs

cd sample_outputs

Rscript ../RSCanner_dotbracket_shannon.R ../sample_test_files/HCV_Jc1.dot ../sample_test_files/HCV_shannon.txt

```

