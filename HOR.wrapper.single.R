#!/usr/bin/env Rscript

#to compile C code on the HPC cluster:

# module load gcc/9
# gcc --version
# gcc -o HOR.V3.3 main.c

# module load mafft/7.309

# module load R/4.0.3



# to be run on the repeats output from TRASH. The wrapper will reformat the repeat files and align using mafft so the C script will work properly. Later, the outputs are read back and HOR summary is produced:
# - repeat weighted consensus
# - HOR repetitiveness within species
# - 

#only repeats that have a class assigned are being used

.libPaths("/home/pw457/Rpackages")

library(stringr)
library(base)
library(msa)#
library(Biostrings)
library(seqinr)#
library(doParallel)#




###########
#Functions#
###########

revCompString = function(DNAstr) {return(toupper(toString(reverseComplement(DNAString(DNAstr))))) }

temp.folder =  "/home/pw457/HOR/singleChromosomes" # "/home/pw457/HORs/temp"

cen180.files.directory = "/home/pw457/HOR" # "/home/pw457/RepeatIdentification/Save/aThaSave"

cen180.files.type = "all.repeats.from" # Possible types: "all.repeats.from", 

HOR.building.strategy = "individual" # Possible types: "individual", 

cutoff = 3
threshold = 5

split.no = 1

task_id = Sys.getenv("SLURM_ARRAY_TASK_ID")
print(paste("slurm array ID for this job is ", task_id))

if(cen180.files.type == "all.repeats.from")
{
  setwd(temp.folder)
  
  repeats.files = list.files(path = cen180.files.directory, pattern = "all.repeats.from", recursive = FALSE) 
  
  if(HOR.building.strategy == "individual")
  {
    pairs.to.compare = data.frame(n = vector(mode = "numeric", length = 0), m = vector(mode = "numeric", length = 0))
    l = 1
    for(n in 1 : 5)
    {
      for(m in n : 5)
      {
        pairs.to.compare = rbind(pairs.to.compare, c(n,m))
        l = l + 1
      }
    }
    names(pairs.to.compare) = c("m", "n")
    
    #foreach(i = 1 : length(repeats.files)) %dopar% {
    #for(i in 1 : length(repeats.files)){
      i = as.numeric(task_id)
      
      repeats = read.csv(paste(cen180.files.directory, repeats.files[i], sep = "/"), stringsAsFactors = FALSE)
      repeats$class[is.na(repeats$class)] = ""
      repeats = repeats[repeats$class == "cen180",]
      repeats$strand[repeats$strand == "+"] = "1"
      repeats$strand[repeats$strand == "-"] = "2"
    
      for(j in 1 : 15){
      print(paste("working on ", j, " pair from ", repeats.files[i], sep = ""))
      
      if(pairs.to.compare$m[j] == pairs.to.compare$n[j])
      {
        repeatsT = repeats[repeats$chromosome == unique(repeats$chromosome)[pairs.to.compare$m[j]],]
        
        file.to.align = paste("to_align_", unique(repeats$chromosome)[pairs.to.compare$m[j]], "_", repeats.files[i], ".fasta", sep = "")
        
        write.fasta(sequences = as.list(repeatsT$seq), names = as.list(paste(1:length(repeatsT$strand), "_D", repeatsT$strand,  sep = "")), file.out = file.to.align, open = "w")
        
        output.alignment = paste("aligned_", unique(repeats$chromosome)[pairs.to.compare$m[j]], "_", repeats.files[i], ".fasta", sep = "")
        
        system(paste("mafft --quiet --inputorder --kimura 1 --retree 1 ", file.to.align, " > ", output.alignment, sep = ""), intern = TRUE)
        
        system(paste("/home/pw457/HOR/HOR.V3.3", " ", temp.folder, "/ ", output.alignment, " ", 1, " ", threshold, " ", cutoff, " ", 1, sep = ""), intern = FALSE)
        
        hors = read.csv(file = paste(temp.folder, "/", "HORs_method_1_", output.alignment, "_t_", threshold, "_c_", cutoff, ".csv", sep = ""), header = TRUE)
        
        if(nrow(hors) > 1)
        {
          hors = as.data.frame(hors)
          
          hors$start.A.bp = repeatsT$start[hors$start_A]
          hors$start.B.bp = repeatsT$start[hors$start_B]
          hors$end.A.bp = repeatsT$end[hors$end_A]
          hors$end.B.bp = repeatsT$end[hors$end_B] 
          hors$chrA = repeatsT$chromosome[hors$start_A]  
          hors$chrB = repeatsT$chromosome[hors$start_B]
          hors$block.size.in.units = hors$end_A - hors$start_A + 1
          hors$block.A.size.bp = hors$end.A.bp - hors$start.A.bp + 1
          hors$block.B.size.bp = hors$end.B.bp - hors$start.B.bp + 1
          for(p in nrow(hors) : 1)
          {
            block1size = abs(hors$end.A.bp[p] - hors$start.A.bp[p])
            block2size = abs(hors$end.B.bp[p] - hors$start.B.bp[p])
            if(abs(block2size - block1size) > 1000)
            {
              hors = hors[-p,]
              p = p - 1
            }
          }
          
          png(filename = paste(temp.folder, "/HOR_plot_", unique(repeats$chromosome)[pairs.to.compare$n[j]], "_", repeats.files[i], ".png", sep = ""), width = 5000, height = 5000, pointsize = 25)
          plot(xlab = unique(repeats$chromosome)[pairs.to.compare$m[j]], ylab = unique(repeats$chromosome)[pairs.to.compare$n[j]], 
               x = hors$start.A.bp, pch = 19, cex = 0.1, 
               y = hors$start.B.bp, 
               xlim = c(min(hors$start.A.bp), max(hors$start.A.bp)), 
               ylim = c(min(hors$start.B.bp), max(hors$start.B.bp)))
          dev.off()
        }
   
        
        remove(hors, repeatsT)
        gc()
      } else
      {
        
      }
    }  
  } else
  {
    
  }
}
print("done")



