# TRASH: Tandem Repeat Annotation and Structural Hierarchy
A package to Identify identify and extract tandem repeats in genome sequences and investigate their higher order structures.

### Overview:
TRASH performs local  **kmer** counting to find regions that are repetitive. Windows are scored based on the proportion of repeated k-mers. Those windows which score above the  **threshold** are considered to contain repeats. The periodicity of tandem repeats is then established and used to find a consensus sequence and map the repeat unites. Periodicity of repeats is established and used to find a consensus sequence and map the repeats. Higher order repeats can be identified as pairs of repeat blocks which are highly similar to each other.

## Requirements:

1. Linux OS (for Windows see below)
2. R-4.1.3 or newer (any R.4+ version should work)

## Installation:

### R
R can be downloaded from https://cloud.r-project.org/ using the instructions provided. Alternatively, a conda enviroment can be set up and activated with:
```
conda create -n name -c conda-forge r-base=4.1.3 zlib
conda activate name
```
### TRASH installer:
Download TRASH and run TRASH_install.sh. This will allow to control whether R packages will be downloaded to a system-default directory or TRASH directory (as in the pre-installed version). Downloaded libraries will be of a specific version, which might cause problems if other versions are already installed.

Adding --def flag to the TRASH_install.sh command will use the default R library path to install new packages. However, note that this will require the use of –def flag each time TRASH is run! This option is suitabel for users that already have multiple R libraries installed and wish to ensure the installation takes as little space as possible by avoiding redundant packages.
```
git clone https://github.com/vlothec/TRASH
cd TRASH
chmod +x TRASH_install.sh
TRASH_install.sh
```
Running TRASH_install.sh can be used to check the installation after it's unpacked and installed
### (alternative) All-in-one package:
Download and unpack TRASH_v1.1.tar.gz which will contain the pre-installed image with all required dependencies:
```
wget https://github.com/vlothec/TRASH/raw/main/TRASH.v1.1.tar.gz
tar -xzvf TRASH.v1.1.tar.gz
```

## Run
TRASH requires at least one fasta file as an input (with ".fa", ".fna" or ".fasta" file extensions). Multiple files can be provided as separate arguments, or by merging sequences into a single fasta file. There is no limit on the number of sequences provided.
### Simple run:
```
TRASH_run.sh assembly.fa --o output.path
```
This will run TRASH with default settings in the output.path directory. 
### Example run
The **/example_run** folder contains a test fasta sequence and results of a TRASH run on this sequence. The sequence is an extraction from chromosome 10 of the CHM13 human genome (coordinates 39,050,443:39,150,442 bp). The test sequence includes several alpha satellite tandem repeat arrays.

## Additional options:

```
--def 			# use the default R packages path.
--rmtemp 		# remove the "*_out" directory after run completion.
--horclass name		# set the name of the repeat family that should be used for HOR calculations, required for the HOR module to be activated.
--limrepno x		# limit alignment sizes (in bp of total sequence) used during the run to calculate consensus, samples repeats to avoid large alignment operations. 78000 by default
--horonly x		# skip the repeat identification if performed earlier and only calculate HORs, needs to be used together with -horclass flag.
--minhor x		# HORs shorter than this value will be discarded, 3 by default.
--maxdiv x		# pair of repeats with a divergence score higher than this value will not be considered as a potential HOR, 5 by default.
--maxchr x  		# the total number of sequences that should be analysed. Useful when the assembly contains a large number of contigs. Sequences are chosen based on their size.
--k x			# kmer size, 10 by default. Decrease if more degraded arrays should be identified, increase for extra stringency (range of 8-16 recommended).
--t x 			# threshold score to choose repetitive windows, 5 by default. Varying this will work similar to k-mer size changes.
--win x 			# window size to use for initial count of repeat content, 1000 by default. Identified repeats will not be larger than this value.
--m x 			# max repeat size to be identified, hard capped by -win setting.
--freg x 		# regions smaller than this will be filtered out at initial steps (some might remain if they come from splitting of a larger region).
--frep x 		# repeats shorter than this will be filtered out, 4 by default.
--o path			# output path where repeats will be saved and temporary directories created.
--seqt path 		# path to the file with repeat family templates, the file needs to be formatted as described below.
--par x 			# max number of cores used for multithreading, defaults to 1. If set as 0, TRASH will try to register as many cores as there are sequences, or maximum available, whichever is smaller.
--randomseed x		# set a random seed for reproducibility of the repeat identification, seed from the previous run can be found in “TRASH_YYYYMMDDHHMMSS.out”.
--simpleplot		# output a plot with repeat coordinates and their sizes for each sequence (additionally to the circos plot).

```


## Multithreading
The script will utilize a maximum of 1 core per fasta sequence (not per file) if available. By default it will use 1 core, which can be controlled with -par flag. 


## Higher Order Repeat (HOR) analysis
TRASH is able to calculate HORs defined as multi-monomer repeat duplications. It does not try to create a 1-dimentional description of repeat monomers, but uses a 2-dimentional matrix of identity between repeats to find instances of consecutive rows of high similarity. -minhor and -maxdiv control how many repeats constitute a HOR and what is the maximum divergence score between repeats for them to be part of a HOR.


## Sequence templates
An additional .csv file can be provided for the run that contains information on predicted repeat families (here called ’class‘). TRASH will check against these templates and if it finds a match, repeats of the same family will be tagged with the provided name. The csv file consists of 3 columns with the names of "name", "length" and "seq". An example file for Arabidopsis thaliana CEN178 would look like:
```
name,length,seq
CEN178,178,AGTATAAGAACTTAAACCGCAACCGATCTTAAAAGCCTAAGTAGTGTTTCCTTGTTAGAAGACACAAAGCCAAAGACTCATATGGACTTTGGCTACACCATGAAAGCTTTGAGAAGCAAGAAGAAGGTTGGTTAGTGTTTTGGAGTCGAATATGACTTGATGTCATGTGTATGATTG
```

## Output

**"all.repeats.from.assembly.fa.csv"**
A table with the tandem repeats identified, their start and end positions, class (when applicable, if not assigned this will be "NA"), sequence (on the positive strand), and strand information (when repeats are assigned to a class they will be identified according to the provided template, thus possibly placing them on the negative strand. In this case, an additional column will contain sequence information on the negative strand).

**"TRASH_assembly.fa.gff"**
GFF file (https://www.ensembl.org/info/website/upload/gff.html) containing the same information as the .csv repeats file but in a format that can be widely used for genome annotation.

**"plots/assembly_circos.pdf"**
A circos plot graphically showing information contained in the "RepetitiveRegions_assembly.fa.csv" file plotted against the analysed sequence.

**"Summary.of.repetitive.regions_assembly.fa.csv"**
Table with regions (arrays) containing tandem repeats, that includes information on their consensus sequence, period size and class (family) when sequence templates are provided

**temp.all.repeats.from.assembly.fa.csv**
A temporary file with repeats, which can be safely removed once the main file is created. This file can be used for troubleshooting.

Optional: **"HOR/HOR_plot_assembly.fasta_chrName_class.png"**
When HORs are calculated, these files are dot plots showing the start locations of the HOR blocks identified.

Optional: **"HORs_assembly.fasta_chrName_class.png.csv"**
A table listing the HORs identified. Each row contains a pair of HOR blocks, each with their start and end coordinate, the number of divergent positions between them (“total variant”) and the direction (1 means the blocks are in the same orientation, i.e. "head to tail", while 2 means they are on opposite strands, "head to head").

**/assembly_out** Specifies the directory with temporary files used during the run that can be removed
**/plots** Specifies the directory containing circos plots, edit.distance plots and HOR.score plots 
Optional: **/HOR** directory with Higher Order Repeat files


## Windows
Windows functionality has not been fully tested (HOR module)
```
git clone https://github.com/vlothec/TRASH
```
extract
Identify Rscript.exe directory
navigate to TRASH\src directory
install TRASH with:
```
[R installation directory]\Rscript.exe TRASH_install.R
```
run TRASH with:
```
[R installation directory]\Rscript.exe TRASH_run.R [run commands]
```
libs.zip are pre-installed Windows libraries that can be downloaded separately