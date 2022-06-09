# RSCanner

RSCanner is an R package for rapid assessment and visualization of RNA structural content that is particularly useful for long RNAs. The script has been written as an Rmd file, and command line functionality has been provided via two separate R scripts, to be used depending on the type of input given.

## file inputs and format

RSCanner takes two types of input containing secondary structural information for the target RNA: 1 - RNA secondary structure in either dot-bracket notation (FASTA file) or CT (connectivity table); 2 - positional Shannon entropy for each nucleotide in the RNA, calculated from base pairing probabilities. Both inputs can be obtained from secondary structure prediction programs (e.g., SuperFold (Smola et al. 2015), Fold from RNAstructure (Reuter and Matthews, 2010), RNAfold from ViennaRNA package (Lorenz et al. 2011)).

CT file should be formatted in the following way: first line is a header "information" line, then six columns of CT information.

dotbracket FASTA file should be formatted in the following way: the file must have three lines - title line, sequence line, dotbracket line, and no newline characters except for the ones separating these three lines.

shannon entropy TEXT file should be formatted in the following way: two columns in a tab delimited file, col1 = index, col2 = shannon entropy values with no header

## usage: RSCanner_CT_shannon.R
To be used when the secondary structure input is a CT file.
For usage information and detailed input/output information, run the following in your terminal:

```
 Rscript RSCanner_CT_shannon.R
```


## usage: RSCanner_dotbracket_shannon.R
To be used when the secondary structure input is in the dot bracket format (FASTA file).
For usage information and detailed input/output information, run the following in your terminal:

```
 Rscript RSCanner_dotbracket_shannon.R
```

## Example
Here is a fully worked example using the sample CT, dotbracket, and shannon text files included in this repository. Simply download the entire repository as-is, and run the following commands from within this repo directory.

### using RSCanner_CT_shannon.R 
```
mkdir sample_outputs

cd sample_outputs

Rscript ../RSCanner_CT_shannon.R ../sample_test_files/HCV_Jc1.ct ../sample_test_files/HCV_shannon.txt

```

### using RSCanner_dotbracket_shannon.R
```
mkdir sample_outputs

cd sample_outputs

Rscript ../RSCanner_dotbracket_shannon.R ../sample_test_files/HCV_Jc1.dot ../sample_test_files/HCV_shannon.txt

```

