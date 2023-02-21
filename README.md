# TRASH: Tandem Repeat Annotation and Structural Hierarchy
Identify and extract tandem repeats and investigate their higher order structure 

### Overview:
Local **kmer** counting finds regions that are repetitive. Windows are scored based on the proportion of repeated kmers to their size. Those which score above the **threshold** are considered to contain repeats. Periodicity of repeats is established and used to find a consensus sequence and map the repeats.

## Requirements:

1. Linux OS
2. R-4.1.3 or newer

## Installation:

### R
R can be downloaded from https://cloud.r-project.org/ using the instructions provided, alternatively a conda enviroment can be set up and activated with:
```
conda create -n name -c conda-forge r-base=4.1.3 zlib
conda activate name
```
### Quick and easy:
Download and unpack TRASH_v1.0.pck.tar.gz which will contain pre-installed image with all required dependancies:
```
wget https://github.com/vlothec/TRASH/raw/main/TRASH_v1.1.tar.gz
tar -xzvf TRASH_v1.0.pck.tar.gz
```
### A bit longer:
Download and unpack TRASH.v1.0.zip and run TRASH_install.sh. This will allow to control whether R packages will be downloaded to a system-default directory or TRASH directory (as in the pre-installed version). Downloaded libraries will be of specific version, which might cause problems if other versions are already installed.

Adding --def flag to the TRASH_install.sh command will use the default R library path to install new packages. This will force the use of --def flag each time TRASH is run!
```
wget --no-check-certificate --content-disposition https://github.com/vlothec/TRASH/TRASH_v1.0.pck.tar.gz
wget https://github.com/vlothec/TRASH/raw/main/TRASH_v1.0.pck.tar.gz
gunzip TRASH.v1.0.zip
chmod +x TRASH_install.sh
TRASH_install.sh
```

## Run
TRASH requires at least one fasta file as an input (with ".fa", ".fna" or ".fasta" extensions). Multiple files can be provided as separate arguments or by merging sequences into one fasta file. There is no limit on the amount of sequences provided.
### Simple run:
```
TRASH_run.sh assembly.fa
```
This will generate 3 files: "RepetitiveRegions_assembly.fa.csv", "Repeats_assembly.fa.csv" and "Repeats_assembly.fa.gff" in the directory from which the command was run. Additionally, 3 directories will be created: "plots_assembly.fa" with circos plots and "assembly.fa_out" with temporary files that can be removed.

## Additional options:

```
--def 			# use the default R packages path.
--rmtemp 		# remove the "*_out" directory after run completion.
--horclass name		# set the name of the repeat family that should be used for HOR calculations, required for the HOR module to be activated.
--limrepno x		# limit alignment sizes (in bp of total sequence) used during the run to calculate consensus, samples repeats to avoid large alignment operations. 78000 by default
--horonly x		# skip the repeat identification if performed earlier and only calculate HORs, needs to be used togehter with -horclass flag.
--minhor x		# HORs shorter than this value will be discarded, 3 by default.
--maxdiv x		# pair of repeats with divergence score higher than this value will not be considered as a potential HOR, 5 by default.
--maxchr x  		# total number of sequences that should be analysed. Usefull when assembly contains large number of contigs. Sequences are chosen based on their size.
--k x			# kmer size, 10 by default. Decrease if more degraded arrays should be identified, increase for extra stringency (range of 8-16 recommended).
--t x 			# threshold score to choose repetitive windows, 5 by default. Change will work similar to the kmer size changes.
--win x 			# window size to use for initial count of repeat content, 1000 by default. Identified repeats will not be bigger than this value.
--m x 			# max repeat size to be identified, hard capped by -win setting.
--freg x 		# regions smaller than this will be filtered out at initial steps (some might remain if they come from splitting of a larger region).
--frep x 		# repeats shorter than this will be filtered out, 4 by default.
--o path			# output path where repeats will be saved and temporary directories created.
--seqt path 		# path to the file with repeat family templates, the file needs to be formatted as described below.
--par x 			# max number of cores used for multithreading, defaults to 1. If set as 0, TRASH will try to register as many cores as there are sequences, or maximum available, whatever is smaller.
--randomseed x		# set a random seed for reproducibility of the repeat identification, seed from the previous run can be found in "TRASH_YYYYMMDDHHMMSS.out" 
--simpleplot		# output a plot with repeat coordinates and their sizes for each sequence (additionally to the circos plot)
```


## Multithreading
The script will utilize a maximum of 1 core per fasta sequence (not per file) if available. By default it will use 1 core, which can be controlled with -par flag. 


## Higher Order Repeat analysis
TRASH is able to calculate HORs defined as multi-monomer repeat duplications. It does not try to create a 1-dimentional description of repeat monomers, but uses a 2-dimentional matrix of identity between repeats to find instances of consecutive rows of high similarity. -minhor and -maxdiv control how many repeats constitute a HOR and what is the maximum divergence score between repeats to be part of a HOR.


## Sequence templates
An additional .csv file can be provided for the run that contains information on predicted repeat families (here called "class"). TRASH will check against it and if it finds a match, repeats of the same family will be tagged with the provided name. It consists of 3 columns with names of "name", "length" and "seq". An example file for Arabidopsis thaliana CEN180 would look like:
```
name,length,seq
CEN180,178,AGTATAAGAACTTAAACCGCAACCGATCTTAAAAGCCTAAGTAGTGTTTCCTTGTTAGAAGACACAAAGCCAAAGACTCATATGGACTTTGGCTACACCATGAAAGCTTTGAGAAGCAAGAAGAAGGTTGGTTAGTGTTTTGGAGTCGAATATGACTTGATGTCATGTGTATGATTG
```

## Output

**"RepetitiveRegions_assembly.fa.csv"**
Table with regions (arrays) containing tandem repeats with information on their consensus sequence, period size and class (family) when sequence templates were provided.

**"Repeats_assembly.fa.csv"**
Table with repeats identified, their start and end positions, class (when applicable, if not assigned it will be "NA"), sequence (on the positive strand) and strand information (when repeats are assigned to a class they will be identified according to the provided template, thus possibly placing them on the negative strand, in this case an additional column will contain sequence information on the negative strand).

**"Repeats_assembly.fa.gff"**
GFF file (https://www.ensembl.org/info/website/upload/gff.html) containing the same information as the .csv repeats file but in a format that can be widely used for genome annotation.

**"plots/assembly_circos.pdf"**
Circos plot showing information contained in the "RepetitiveRegions_assembly.fa.csv" file.

Optional: **"HOR/HOR_plot_assembly.fasta_chrName_class.png"**
When HORs were calculated, these are dot plots showing the start locations of HOR blocks identified.

Optional: **"HORs_assembly.fasta_chrName_class.png.csv"**
Table with HORs identified. Each row is a pair of HOR blocks, each with their start and end coordinate, number of divergent positions between them ("total variant") and direction (1 means the blocks are in the same orientation, i.e. "head to tail", while 2 means they are on opposite strands, "head to head").
