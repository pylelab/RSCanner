# RSCanner

RSCanner is an R package for rapid assessment and visualization of RNA structural content that is particularly useful for long RNAs. The script has been written as an Rmd file, and command line functionality has been provided via two separate R scripts, to be used depending on the type of input given.

## Inputs and algorithm description

Input files: 1 - RNA secondary structure in either dot-bracket notation (FASTA file) or CT (connectivity table); 2 - positional Shannon entropy for each nucleotide in the RNA, calculated from base pairing probabilities.

RSCanner scans the secondary structure provided along with Shannon entropy values in sliding windows covering the entire RNA and calculates the base pair content (BPC) and median entropy values for each window. Then, it computes nucleotide positions with BPC values above a user-defined cutoff (default = 50th percentile) and Shannon entropy values below a user-defined cutoff (default = 50th percentile). These positions are termed structure counts. RSCanner then plots the frequency distribution of structure counts across the RNA as heatmap and histogram plots.  
How are ends handled? Include here.

## Usage: RSCanner_CT_shannon.R
To be used when the secondary structure input is a CT file.
For usage information and detailed input/output information, run the following in your terminal:

```
 Rscript RSCanner_CT_shannon.R
 
 # Terminal output
 
 To run the script, do 
 Rscript RSCanner_CT_shannon.R <input_structure.ct>  <shannon.txt>

 Input:
    input_structure.ct             - CT format file; six columns after one header line
    shannon.txt               - two columns in a tab delimited file, col1 = index, col2 = Shannon Entropy values with no header 

 Output:
    bpcplot.tiff                             - output line plot figure
    smoothed_shannonplot.tiff                - output line plot figure
    structure_counts_histogram.tiff          - output histogram figure
    structure_counts_heatmap.tiff            - output heatmap figure
    
    bpc_data.csv                             - output tab delimited file, col1=index, col2=necleotide number, col3=BPC
    smoothed_Shannon_data.csv                - output tab delimited file, col1=index, col2=necleotide number, col3=smoothed Shannon Entropy
    structure_counts.csv                     - output tab delimited file, col1=index, col2=structrue counts, col3=BPC, col4=smoothed Shannon Entropy
    ordered_structure_table.csv              - output tab delimited file, col1=index, col2=Bin Number, col3=%structure counts, col4=Bin start(nt), col5=Bin end(nt) 


```

## Usage: RSCanner_dotbracket_shannon.R
To be used when the secondary structure input is in the dot bracket format (FASTA file).
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

