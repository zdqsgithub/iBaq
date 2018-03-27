Intensity-based Absolute Quantification (iBaq) r Script
DQZ, Graham Lab, USC

Future update:

1.Merge in-silico digested iBaq normalization factors proteome library script

2.Fisher's combined method to calculate protein level p-value

3.create user-friendly R package


03/26/2018
iBaq v1.1 r script update:

1.For peptides that are mapped to multiple proteins (A; B; C), split these rows into individual rows, one for each protein mapping (A; B; C becomes three rows, one for A, one for B, one for C)

2.Keep track of the “uniqueness” of each peptide → protein mapping.
    -A peptide that maps to only one protein scores 1
    -A peptide that maps to multiple proteins scores 0

3.Using all the peptide → protein mappings, do iBAQ calculation. At the same time, for each protein that emerges from iBAQ, calculate:
    -The number of unique peptides used for the iBAQ value (ie, the sum of the “uniqueness” score)
    -The total number of peptides used for the iBAQ value (ie, the total number of peptides = the number of peptides uniquely mapped to that protein + the number of peptides mapped to that protein and also other proteins)

4.The unique and multiple mapped peptides scores for each iBAQ protein can then be used for later filtering

5.uploaded a testing "test-in.txt" file (500 peptides) and added tutorial instruction in the script comments

6.change variable to meaningful names 

7.The script is simplify by calling the function: pepGroup_to_iBaq(Peptides,sampleStart,sampleEnd,samples,outputFileName)

8.validation by hand calculation




03/23/2018 update:
1. perform iBaq calculation PD 2.1 Minora LFQ peptideGroups.txt

2.iBaq-v1.0 requires in-silico digested iBaq normalization factors library

3."20170822-Human-uniprot-all-reviewed-iBAQ-counter-only.txt" contains Human-uniprot-all-reviewed iBaq normalization factors, calculated by a separate script (will update)

4.pushed existing iBaq v1.0 r script to GitHub repository


