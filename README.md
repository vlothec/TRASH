## RepeatIdentifier.R
The first part of the script uses a fasta file as input and finds tandemly arranged repeats.
## HOR.wrapper.single.R
The second part of the script, UNFINISHED, uses a .csv file with repeats and finds higher order repeats.
## HOR.V3.3.c 
C script required for the HOR part

# TRASH: Tandem Repeat Annotation and Structural Hierarchy
identify and extract tandem repeats and investigate their higher order structure 

Local **kmer** counting finds regions that are repetitive. Windows (1 kbp by default) are scored based on the proportion of repeated kmers identified in relation to their size. **Threshold** is a score above which windows are considered to contain repeats. The periodicity of repeats is established and **MAFFT** is used to iterate over candidate representative sequences to find a consensus sequence. 

## Before running:

1. Make sure you have required R packages installed. **Required**
2. Create a directory to write the ouptuts to. Set the path using the "outputs.directory" variable. **Required**
3. Choose/create a directory containing the fasta formatted sequences to be analysed. Set the path to the directory using the "genomes.directory" variable. **Required**
4. Prepare a .csv file containing repeat templates that can be used to shift the identified repeats such that their start sites align. Set its path using the "sequence.templates" variable. More information related to this is located in the script's comments **Optional**
5. The script automatically tries to use as many CPUs as there are DNA sequences in the fasta files in the provided directory. If less cores should be used, the "set.no.of.cores" variable can be changed. Multithreading is not implemented for Windows OS **Optional**

## Outputs:
For each fasta file input there will be created a separate directory containing:
1. frame4 and frame2 directories that contain alignments. These can take up a lot of space and can be removed.
2. plots directory contains plots of each sequence in the fasta file. The plots are of repeat size (y axis) and position along the analysed target sequence (x axis), plotted for each repeat.
3. All.repeats.from.NAME.csv file that contains all repeats in a .csv format.
4. Summary.of.repetitive.regions.NAME.csv file that contains details of the regions that were analysed.


## Other settings:
TODO

## Misc
The script was tested on Windows and Linux OS

# HOR 
TODO


