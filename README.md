# TRASH: Tandem Repeat Annotation and Structural Hierarchy
Identify and extract tandem repeats and investigate their higher order structure 

### Ooverwiev:
Local **kmer** counting finds regions that are repetitive. Windows (1kb by default) are scored based on proportion of repeated kmers to their size. Windows above the **threshold** are considered to contain repeats. Periodicity of repeats is established and used to find a consensus sequence and map the repeats.

## Requirements:

1. Linux-OX
2. R-4.1.3 or newer

## Installation:

### R
R can be downloaded from https://cloud.r-project.org/ using the instructions provided, alternatively a conda enviroment can be set up and activated with:
```
conda install -n name -c conda-forge r-base=4.1.3
conda activate name
```
### Quick and easy:
Download and unpack TRASH.v1.0.pck.zip which will contain pre-installed image with all required dependancies:
```
wget --no-check-certificate --content-disposition https://github.com/vlothec/TRASH/TRASH.v1.0.pck.zip
gunzip TRASH.v1.0.pck.zip
chmod +x TRASH_run.sh
```
You can add TRASH_run.sh to the PATH directory
```
export PATH=$PATH:/TRASH_dir/TRASH_run.sh
```
### A bit longer:
Download and unpack TRASH.v1.0.zip and run TRASH_install.sh. This will allow to control whether R packages will be downloaded to a system-default directory or TRASH directory (as in the pre-installed version). Downloaded libraries will be of specific version, which might cause problems if other versions are already installed.

Adding --def flag to the TRASH_install.sh command will use of default R library path to install new packages. This will force the use of --def flag each time TRASH is run!
```
wget --no-check-certificate --content-disposition https://github.com/vlothec/TRASH/TRASH.v1.0.pck.zip
gunzip TRASH.v1.0.zip
chmod +x TRASH_install.sh
TRASH_install.sh
```

## Run
TRASH requires at least one fasta file as an input (with ".fa", ".fna" or ".fasta" extensions). Multiple files can be provided as separate arguments or by merging sequences into one fasta file. There is no limit on amount of sequences provided.
### Simple run:
```
TRASH_run.sh assembly.fa
```
This will generate 3 files: "RepetitiveRegions_assembly.fa.csv", "Repeats_assembly.fa.csv" and "Repeats_assembly.fa.gff" in the directory from which the command was run. Additionally, 3 directories will be created: "plots_assembly.fa" with circos plots and "assembly.fa_out" with temporary files that can be removed.

## Additional options:

### Repeat identification settings
```
-def 			#use the default R packages path.
-rmtemp 		#remove the "*_out" directory after the run completion.
-horclass name	#set the name of the repeat family that should be used for HOR calculations, required for the HOR module to be activated.
-limrepno x			#limit alignment sizes (in bp of total sequence) used during the run to calculate consensus, samples repeats to avaoid large alignment operations. 78000 by default
-horonly x		#skip the repeat identification if was performed earlier and only calculate HORs, needs to be used togehter with -horclass flag.
-minhor x		#HORs shorter than this value will be discarded, 3 by default.
-maxdiv x		#pair of repeats with divergence score higher than this value will not be considered as a potential HOR, 5 by default.
-maxchr x  		#total number of sequences that should be analysed. Usefull when assembly contains large number of contigs. Sequences are chosen based in their size.
-k x			#kmer size, 10 by default. Decrease if more degraded arrays should be identified, increase for extra stringency (range of 8-16 recommended).
-t x 			#threshold score to choose repetitive windows, 5 by default. Change will work similar to the kmer size changes.
-win x 			#window size to use for initial count of repeat content, 1000 by default. Identified repeats will not be bigger than this value.
-m x 			#max repeat size to be identified, hard capped by -win setting.
-freg x 		#regions smaller than this will be filtered out at initial steps (some might remain if they come from splitting of a larger region).
-frep x 		#repeats shorter than this will be filtered out, 4 by default.
-o path			#output path where repeats will be saved and temporary directories created.
-seqt path 		#path to the file with repeat family templates, the file needs to be formatted as described below.
-par x 			#max number of cores ussed for multithreading, defaults to 1. If set as 0, TRASH will try to register as many cores as there are sequences, or maximum available, whatever is smaller.
```


## Threading
The script will utilize maximum of 1 core per fasta sequence (not per file) if available. Dy default it will use up to 1 cores, which can be controlled with -par flag. 


## Higher Order Repeat analysis

TRASH is able to calculate HORs defined as multi-monomer repeat duplications. It does not try to create a 1-dimentional description of repeat monomers, but uses a 2-dimentional matrix of identity between repeats to find instances of consecutive rows of high similarity. -minhor and -maxdiv control how many repeats constitute a HOR and what is the maximum divergence score between repeats to be part of a HOR.



