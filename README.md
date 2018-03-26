# iBaq
iBAQ LFQ on PD 2.1 peptideGroups.txt. A Graham Lab USC project. 

03/23/2018
Importing existing iBaq r script
iBaq-v2 reuires in-silico digested iBaq normalization factors proteome library
can perform iBaq calculation PD 2.1 Minora LFQ peptideGroups.txt export 
issue: protein assemble from peptides

Taskes (priority from high to low):

Update iBAQ script:

1.For peptides that are mapped to multiple proteins (A; B; C), split these rows into individual rows, one for each protein mapping (A; B; C becomes three rows, one for A, one for B, one for C)

2.Keep track of the “uniqueness” of each peptide → protein mapping.
    -A peptide that maps to only one protein scores a 1
    -A peptide that maps to multiple proteins scores a 0

3.Using all the peptide → protein mappings, do iBAQ calculation. At the same time, for each protein that emerges from iBAQ, calculate:
    -The number of unique peptides used for the iBAQ value (ie, the sum of the “uniqueness” score)
    -The total number of peptides used for the iBAQ value (ie, the total number of peptides = the number of peptides uniquely mapped to that protein + the number of peptides mapped to that protein and also other proteins)

4.The unique and multiple mapped peptides scores for each iBAQ protein can then be used for later filtering

5.Merge in-silico digested iBaq normalization factors proteome library script

6.create user-friendly R package



