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
Download and unpack TRASH_v1.2.tar.gz which will contain the pre-installed image with all required dependencies:
```
wget https://github.com/vlothec/TRASH/raw/main/TRASH.v1.2.tar.gz
tar -xzvf TRASH.v1.2.tar.gz
```

## Run
TRASH requires at least one fasta file as an input (with ".fa", ".fna" or ".fasta" file extensions). Multiple files can be provided as separate arguments, or by merging sequences into a single fasta file. There is no limit on the number of sequences provided.
### Simple run:
```
TRASH_run.sh assembly.fa --o output.path
```
This will run TRASH with default settings in the output.path directory. 
Absolute path required at this point. Not specifying path at all will direct the outputs to the curent directory
### Example run
The **/example_run** folder contains a test fasta sequence and results of a TRASH run on this sequence. The sequence is an extraction from chromosome 10 of the CHM13 human genome (coordinates 39,050,443:39,150,442 bp). The test sequence includes several alpha satellite tandem repeat arrays.

## Additional options:

```
--def 			# use the default R packages path.
--rmtemp 		# remove the "*_out" directory after run completion.
--simpleplot		# output a plot with repeat coordinates and their sizes for each sequence (additionally to the circos plot).
--horclass name		# set the name of the repeat family that should be used for HOR calculations, required for the HOR module to be activated.
--limrepno x		# limit alignment sizes (in bp of total sequence) used during the run to calculate consensus, samples repeats to avoid large alignment operations. 78000 by default
--horonly x		# skip the repeat identification if performed earlier and only calculate HORs, needs to be used together with -horclass flag.
--minhor x		# HORs shorter than this value will be discarded, 3 by default.
--maxdiv x		# pair of repeats with a divergence score higher than this value will not be considered as a potential HOR, 5 by default.
--maxchr x  		# the total number of sequences that should be analysed. Useful when the assembly contains a large number of contigs. Sequences are chosen based on their size.
--k x			# kmer size, 10 by default. Decrease if more degraded arrays should be identified, increase for extra stringency (range of 8-16 recommended).
--t x 			# threshold score to choose repetitive windows, 5 by default. Varying this will work similar to k-mer size changes.
--win x 		# window size to use for initial count of repeat content, 1000 by default. Identified repeats will not be larger than this value.
--m x 			# max repeat size to be identified, hard capped by -win setting.
--freg x 		# regions smaller than this will be filtered out at initial steps (some might remain if they come from splitting of a larger region).
--frep x 		# repeats shorter than this will be filtered out, 4 by default.
--o path		# output path where repeats will be saved and temporary directories created. Use an absolute path
--seqt path 		# path to the file with repeat family templates, the file needs to be formatted as described below.
--par x 		# max number of cores used for multithreading, defaults to 1. If set as 0, TRASH will try to register as many cores as there are sequences, or maximum available, whichever is smaller.
--randomseed x		# set a random seed for reproducibility of the repeat identification, seed from the previous run can be found in “TRASH_YYYYMMDDHHMMSS.out”.
--N.max.div x       	# (monomer splitting method) threshold score above which will look for divisions, the lower, the more loose. 100 by default, meaning the method is turned off. Suggested setting when a monomer merging arise is 5.
--max.N.split x         # (monomer splitting method) max number of N divisions, the higher the longer repeats can be split. 12 by default
--smooth.percent x  	# (monomer splitting method) smoothing factor for finding the histogram peaks, the higher the wider. 2 by default

```


## Multithreading
The script will utilize a maximum of 1 core per fasta sequence (not per file) if available. By default it will use 1 core, which can be controlled with --par flag. 

## Higher Order Repeat (HOR) analysis
TRASH is able to calculate HORs defined as multi-monomer repeat duplications. It does not try to create a 1-dimentional description of repeat monomers, but uses a 2-dimentional matrix of identity between repeats to find instances of consecutive rows of high similarity. --minhor and --maxdiv control how many repeats constitute a HOR and what is the maximum divergence score between repeats for them to be part of a HOR.

TRASH will also calculate the sum of the lengths of all HORs each repeat is a part of and report it in the repeats table under "repetitiveness" column. 

## Sequence templates
An additional .csv file can be provided for the run that contains information on predicted repeat families (here called ’class‘). TRASH will check against these templates and if it finds a match, repeats of the same family will be tagged with the provided name. The csv file consists of 3 columns with the names of "name", "length" and "seq". An example file for Arabidopsis thaliana CEN178 would look like:
```
name,length,seq
CEN178,178,AGTATAAGAACTTAAACCGCAACCGATCTTAAAAGCCTAAGTAGTGTTTCCTTGTTAGAAGACACAAAGCCAAAGACTCATATGGACTTTGGCTACACCATGAAAGCTTTGAGAAGCAAGAAGAAGGTTGGTTAGTGTTTTGGAGTCGAATATGACTTGATGTCATGTGTATGATTG
```
If sequence templates are provided, TRASH is going to align all repeats from each class within each sequence, generate consensus and calculate the Levenshtein distance between each repeat and the consensus as a divergence metric.

## Output

**"all.repeats.from.assembly.fa.csv"**
A table of the tandem repeat monomers identified, their start and end positions, class (when applicable if templates have been used, if no templates were provided it will be assigned ‘NA’), sequence (on the positive strand), and strand information. When repeats are assigned to a class they will be identified according to the provided template, thus possibly placing them on the negative strand, an additional column will contain sequence information on the negative strand. Classified repeats also have their edit distance from the sequence-wide consensus calculated. If the HOR identification module was used, a repetitiveness score which is a sum of all HOR lengths (in monomers) that a repeat is a part of is also reported here. 

**"TRASH_assembly.fa.gff"**
A GFF file containing the same information as the ‘all.repeats.from.assembly.fa.csv’ file, but in a format that can be widely used for genome annotation analysis and visualisation.

**"plots/assembly_circos.pdf"**
Circos plot showing the identified repeat monomers using coordinates in the ‘Summary.of.repetitive.regions_assembly.fa.csv’ file plotted against the fasta sequences analysed.

**"Summary.of.repetitive.regions_assembly.fa.csv"**
A table with regions (arrays) containing tandem repeats, including information on their consensus sequence, monomer size and class (family) when sequence templates are provided.

**temp.all.repeats.from.assembly.fa.csv**
A temporary file containing repeat annotation that can be safely removed once the main file ‘all.repeats.from.assembly.fa.csv’ is created. 

Optional: **"HOR/HOR_plot_assembly.fasta_chrName_class.png"**
When HORs are identified, these files provide dot plots showing the start locations of HOR blocks. 

Optional: **"HORs_assembly.fasta_chrName_class.png.csv"**
A table of identified HORs. Each row reports a pair of HOR blocks, with their start and end coordinates, the number of divergent positions between them (‘total variant’) and direction (1 means the blocks are in the same orientation, i.e. ‘head to tail’, while 2 means they are on opposite strands, ‘head to head’).

**/assembly_out** Specifies the directory with temporary files used during the run that can be removed
**/plots** Specifies the directory containing circos plots, edit.distance plots and HOR.score plots 
Optional: **/HOR** directory with Higher Order Repeat files

## Monomer splitting method
In some cases, the signal from a multiplication of the base repeat might be stronger than the one from the base repeat itself, resulting in identification of repeats in multimers. To address it, TRASH divides the most frequent k-mer N by a range of integers (2 to 12 by default, with the upper limit controlled using the ‘--max.N.split’ flag) and checks whether peaks exist at these new k-mer distances. For each integer d (2 to 12 by default), TRASH will sum k-mers found surrounding the N/d distance, and take the highest possible d that is above a percentage threshold set by the user. This threshold, controlled with the ‘--N.max.div’ flag, is set at 100 by default, meaning the method is normally not functional. When considering composite numbers (4, 6, 8, 9 etc), TRASH will also consider the number of k-mers around distance values that correspond to division of N by all the positive divisors (other than 1 and itself).

## Windows
Windows functionality has not been fully tested (especially the HOR module), but most call-outs and installation have been adjusted for Windows use.

Installation:
```
git clone https://github.com/vlothec/TRASH
```
1. Extract the archive
2. Identify Rscript.exe directory
3. navigate to TRASH\src directory
4. install TRASH with:
```
[R installation directory]\Rscript.exe TRASH_install.R
```
run TRASH with:
```
[R installation directory]\Rscript.exe TRASH_run.R [run commands]
```
win_libs.zip are pre-installed Windows libraries that can be unpacked directly into /libs directory

## Cite

Piotr Wlodzimierz and others, TRASH: Tandem Repeat Annotation and Structural Hierarchy, 
Bioinformatics, Volume 39, Issue 5, May 2023, btad308, 
https://doi.org/10.1093/bioinformatics/btad308

