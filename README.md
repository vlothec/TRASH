# TRASH: Tandem Repeat Annotation and Structural Hierarchy
R script to identify and extract tandem repeats

Local **kmer** counting finds regions that are repetitive. Windows (1kb by default) are scored based on proportion of repeated kmers to their size. **Threshold** is a score above which windows are considered to contain repeats. Periodicity of repeats is established and **MAFFT** is used to iterate over candidate representative sequences to find a consensus. 

## Before running:

1. Make sure you have required R packages installed. **Required**
2. Create a directory to write the ouptuts to, set the path as "outputs.directory" variable. **Required**
3. Choose/create a directory with fasta formatted sequences to be analysed, set the path (to the directory) as "genomes.directory" variable. **Required**
4. Prepare a csv file with repeat templates that can be used to shift the identified ones so that their start sites align. Set its path as "sequence.templates" variable. More on that in the script's comments **Optional**
5. The script automatically tries to use as many CPUs as there are DNA sequences in the fasta files in the provided directory. If less cores should be used, "set.no.of.cores" variable can be changed. Multithreading not implemented for Windows OS **Optional**

## Outputs:
For each fasta file input there will be created a separate directory with:
1. frame4 and frame2 directories that contain alignments. Can take up a lot of space and can be removed
2. plots directory contains plots of each sequence in the fasta file with size (y axis) and position (x axis) plotted for each repeat
3. All.repeats.from.NAME.csv file that contains all repeats in a csv format
4. Summary.of.repetitive.regions.NAME.csv file that contains details of the regions that were being analysed


## Other settings:
TODO

## Misc
The script was tested on Windows and Linux OS


