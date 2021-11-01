#!/usr/bin/env Rscript
library(stringr)
library(base)
library(msa)
library(Biostrings)
library(seqinr)
library(doParallel)


##################
#General settings#
##################

#Multithreading#
#it's based on amount of .fasta files provided as an input, to use more CPUs when running a single genome, it can be split into chromosomes
#if left at 0, R will try to use the same amount or cores as the amount of sequence files in the input folder
set.no.of.cores = 0 

#Run settings#
set.kmer = 12 # kmer size used for initial identification of repetitive regions
set.threshold = 10 # window repetitiveness score (0-100) threshold
set.max.repeat.size = 500 # max size of repeats to be identified
outputs.directory = "" # output folder. example "~/Arabidopsis/output"
genomes.directory = "" # folder with .fasta inputs. example "~/Arabidopsis/genomes"
sequence.templates = NA # path to a csv file used to match identified repeats to templates so they are globally in the same frame (same start position), set as NA to skip
##example: sequence.templates = "~/sequence.template.csv"
##format: 
##seq,name,length,group
##AGTATA,CEN180,180,ath


test.sequence = FALSE # keep FALSE, otherwise outputs and genomes directory paths will be overwritten
###########
#Functions#
###########

revCompString = function(DNAstr) {return(toupper(toString(reverseComplement(DNAString(DNAstr))))) }

extract.kmers = function(fasta.sequence, kmer.size)
{
  fasta.sequence = strsplit(fasta.sequence, split = "")[[1]]
  size = length(fasta.sequence)
  fasta.sequence = c(fasta.sequence, fasta.sequence)
  kmers = vector(mode = "character", length = size)
  for(i in 1 : size)
  {
    kmers[i] = paste(fasta.sequence[i : (i + kmer.size - 1)], collapse = "")
  }
  return(kmers)
}

kmer.compare = function(repA, repB, kmer.size = 8)
{
  if(is.na(repA))
  {
    print("repA not a character!")
    return(0)
  }
  if(is.na(repB))
  {
    print("repB not a character!")
    return(0)
  }
  repAkmers = extract.kmers(repA, kmer.size)
  repBkmers = extract.kmers(repB, kmer.size)
  common = ( length(which(repBkmers %in% repAkmers)) + length(which(repAkmers %in% repBkmers)) ) / 2
  if(is.na(common))
  {
    return(0)
  }
  d = 100 * common/ave(nchar(repA), nchar(repB))
  if(is.na(d))
  {
    return(0)
  }
  return(d)
}

RepeatIdentifier = function(DNA.sequence = "", assemblyName = "", fasta.name = "", 
                            kmer = 12, window = 1000, threshold = 10, mask.small.regions = 1500, mask.small.repeats = 4,
                            max.repeat.size = 500,
                            tests = 5, temp.folder = "",
                            sequence.template)
{
  print(paste("starting the main function", assemblyName, fasta.name, sep = " "))
  regions.data.frame = data.frame(index = integer(), fasta.name = character(), start = integer(), end = integer(), average.score = double(), most.freq.value.N = integer())
  
  index = 0
  
  seq.length = nchar(DNA.sequence)
  
  scores = vector(mode = "numeric", length = seq.length %/% window + 1)
  
  for(ii in 1 : (length(scores) - 1))  # for every window of size window but not the last one
  {
    if(ii%%window == 0)
    {
      print(paste("Checking windows of sequence ", fasta.name, ", window ", ii, "/", length(scores), sep = ""))
    }
    window.start = (ii - 1) * window
    check.window = str_sub(DNA.sequence, window.start, (window.start + window - 1))
    for(iii in 1 : (window - kmer))
    {
      check.string = str_sub(string = check.window, start = iii, end = (iii + kmer - 1))
      count = 1 * ((str_count(string = check.string, pattern = "N") == 0) & str_count(string = check.window, pattern = check.string) > 1)
      scores[ii] = scores[ii] + count
    }
  }
  ii = length(scores) # for the last window
  print(paste("sequence ", fasta.name, ", window ", ii, "/", length(scores), sep = ""))
  window.start = (ii - 1) * window
  check.window = str_sub(DNA.sequence, window.start)
  for(iii in 1 : (nchar(check.window) - kmer))
  {
    check.string = str_sub(string = check.window, start = iii, end = (iii + kmer - 1))
    count = 1 * ((str_count(string = check.string, pattern = "N") == 0) & str_count(string = check.window, pattern = check.string) > 1)
    scores[ii] = scores[ii] + count
  }
  
  scores[1 : (length(scores) - 1)] = scores[1 : (length(scores) - 1)] / window * 100 #change score into a percentage score
  scores[length(scores)] = scores[length(scores)] / nchar(check.window)
  
  ii = 1
  while(ii <= length(scores)) #apply threshold and save repetitive windows from this fasta
  {
    if(scores[ii] >= threshold)
    {
      index = index + 1
      name = fasta.name
      start = (ii - 1) * window
      end = start - 1
      ave.score = 0
      s = 0
      while(scores[ii] >= threshold)
      {
        s = s + 1
        end = end + window
        ave.score = scores[ii] + ave.score 
        ii = ii + 1
      }
      ave.score = ave.score / s
      most.freq.value.N = 0
      temp = data.frame(index, name, start, end, ave.score, most.freq.value.N)
      regions.data.frame = rbind(regions.data.frame, temp)
    }
    ii = ii + 1
  }
  
  regions.data.frame = regions.data.frame[(regions.data.frame$end - regions.data.frame$start) > mask.small.regions,]
  
  if(nrow(regions.data.frame) < 1)
  {
    return(0)
  }
  
  
  print(paste("Step 1 (identify windows) done", sep = ""))
  
  
  for(i in 1 : nrow(regions.data.frame))
  {
    seqB = str_sub(DNA.sequence, regions.data.frame$start[i], regions.data.frame$end[i])
    print(paste("calculating distances on region ", i, "/", nrow(regions.data.frame), " region size: ", nchar(seqB), sep = ""))
    distance = vector(mode = "numeric", length =  nchar(seqB))
    for(ii in 1 : (nchar(seqB) - kmer + 1))
    {
      kmer.pattern = str_sub(seqB, ii, ii + kmer - 1)
      window.string = str_sub(seqB, ii + 1, ii + max.repeat.size)
      distance[ii] = str_locate(string = window.string, pattern = kmer.pattern)[[1]]
    }
    distance = distance[!is.na(distance)]
    distance = distance[distance > 0]
    if(length(distance) > 0)
    {
      a = hist(distance, breaks = max(distance), plot = FALSE)
      regions.data.frame$most.freq.value.N[i] = a$breaks[which.max(a$counts) + 1]
    }
  }
  
  regions.data.frame = regions.data.frame[regions.data.frame$most.freq.value.N >= mask.small.repeats,]
  
  #Test random N-size mers, handle overlaps, choose the best one, extract repeats for mafft alignment
  if(nrow(regions.data.frame) < 1)
  {
    return(0)
  }
  
  consensus.primary = vector(mode = "character", length = nrow(regions.data.frame))
  consensus.count = vector(mode = "numeric", length = nrow(regions.data.frame))
  
  if(temp.folder == "")
  {
    stop("no temp folder")
  }
  setwd(temp.folder)
  
  if(!dir.exists(assemblyName))
  {
    dir.create(assemblyName)
  }
  if(!dir.exists(paste(temp.folder, "/", assemblyName, "/frame2", sep ="")))
  {
    dir.create(paste(temp.folder, "/", assemblyName, "/frame2", sep =""))
  }
  
  for(i in 1 : nrow(regions.data.frame))
  {
    print(paste("testing random N samples from region ", i, "/", nrow(regions.data.frame), " window size: ", regions.data.frame$end[i] - regions.data.frame$start[i], sep = ""))
    seqC = str_sub(DNA.sequence, regions.data.frame$start[i], regions.data.frame$end[i])
    N = regions.data.frame$most.freq.value.N[i]
    if(N == 0)
    {
      N = 1
    }
    random.sequences.start = sample(1 : (nchar(seqC) - N), tests)
    random.sequences = str_sub(seqC, random.sequences.start, (random.sequences.start + N - 1))
    random.sequence.scores = vector(mode = "numeric", length = tests)
    random.sequence.counts = vector(mode = "numeric", length = tests)
    
    maxMis = N %/% 3
    if(maxMis > 100)
    {
      maxMis = 100
    }
    
    #match = countPDict(pdict = random.sequences, subject = seqC, max.mismatch = maxMis, with.indels = TRUE)#check
    
    for(ii in 1 : tests)
    {
      match = matchPattern(pattern = random.sequences[ii], subject = seqC, max.mismatch = maxMis, with.indels = TRUE)
      temp = as.data.frame(match)
      temp$start  = start(match)
      temp$end = end(match)
      temp$strand = "+"
      iii = nrow(temp)
      random.sequence.counts[ii] = 0
      if(iii > 0)
      {
        random.sequence.counts[ii] = nrow(temp) 
      }
      while(iii > 1)
      {
        if(temp$start[iii] + maxMis <  temp$end[iii - 1])
        {
          temp = temp[-iii,]
        }
        iii = iii - 1
      }
      if(nrow(temp) > 1)
      {
        write.csv(temp, file = paste(temp.folder, "/", assemblyName, "/frame2/", fasta.name, "_", i, "_", ii, ".csv", sep = ""), row.names = FALSE)
        random.sequence.scores[ii] = sum(temp$end - temp$start)
      }
    }
    #mafft
    if(random.sequence.scores[which.max(random.sequence.scores)] != 0)
    {
      winner = read.csv(paste(temp.folder, "/", assemblyName, "/frame2/", fasta.name, "_",  i, "_", which.max(random.sequence.scores), ".csv", sep = ""))
      
      if(ncol(winner) == 5)
      {
        write.fasta(file.out = paste(assemblyName, "/frame2/", i, "_", "pre.extract", regions.data.frame$index[i], ".", regions.data.frame$name[i], ".fasta", sep = ""), 
                    names = paste(winner$start, winner$end, sep = "_"), sequences = str_split(winner[,4], pattern = ""), open = "w", as.string = FALSE)
      } else if(ncol(winner) == 4){
        write.fasta(file.out = paste(assemblyName, "/", "/frame2/", i, "_", "pre.extract", regions.data.frame$index[i], ".", regions.data.frame$name[i], ".fasta", sep = ""), 
                    names = paste(winner$start, winner$end, sep = "_"), sequences = str_split(winner$x, pattern = ""), open = "w", as.string = FALSE)
      }
      
      
      input = paste(assemblyName, "/frame2/",  i, "_", "pre.extract", regions.data.frame$index[i], ".", regions.data.frame$name[i], ".fasta", sep = "")
      output = paste(assemblyName, "/frame2/",  i, "_", "pre.extract", regions.data.frame$index[i], ".", regions.data.frame$name[i], ".aligned.fasta", sep = "")
      
      if(Sys.info()['sysname'] == "Linux")
      {
        system(paste("mafft --quiet --retree 2 --inputorder ", input, " > ", output, sep = ""), intern = TRUE)
      } else if(Sys.info()['sysname'] == "Windows")
      {
        shell(paste("mafft --quiet --retree 2 --inputorder ", input, " > ", output, sep = ""), intern = TRUE)
      }
      
      alignment = read.alignment(output, format = "FASTA", forceToLower = FALSE)
      consensus = consensus(alignment, threshold = 0.3)
      consensus = toupper(consensus[consensus != "-"])
      consensus = paste(consensus, collapse = "")
      
      
      consensus.primary[i] = consensus
      consensus.count[i] = nrow(winner)
      
      remove(alignment, input, output, winner)
    } else
    {
      consensus.primary[i] = "none_identified"
      consensus.count[i] = 0
      
    }
  }
  regions.data.frame = cbind(regions.data.frame, consensus.primary)
  regions.data.frame = cbind(regions.data.frame, consensus.count)
  regions.data.frame$fasta.name = fasta.name
  
  print(paste("finished test random Nmers of ", assemblyName, ": ", fasta.name, sep = ""))
  
  write.csv(regions.data.frame, file = paste(assemblyName, "/Pre_shift_regions_", assemblyName, "_", fasta.name, ".csv", sep = ""))
  # 3.5 Shift the primary consensus sequences to match a given sequence, only those regions with N within a given range
  regions.data.frame$class = ""
  if(!is.na(sequence.template)[1])
  {
    for(i in 1 : nrow(regions.data.frame))
    {
      if(regions.data.frame$consensus.primary[i] != "none_identified")
      {
        highest = NULL
        check.seq = regions.data.frame$consensus.primary[i]
        if((nchar(check.seq) > (min(sequence.template$length)-10)) & (nchar(check.seq) < (max(sequence.template$length)+10)))
        {
          scores = NULL
          for(ii in 1 : nrow(sequence.template))
          {
            scores = c(scores, kmer.compare(sequence.template$seq[ii], check.seq))
          }
          highest = which.max(scores)
          
          seq.len = nchar(check.seq)
          #create a df of all possible shifts and a column for their scores
          shifts = data.frame(seq = vector(mode = "character", length = (seq.len * 2)), score = vector(mode = "numeric", length = (seq.len * 2)))
          #split the sequence to make it easier to handle
          seq.to.split.plus = strsplit(check.seq, split = "")[[1]]
          seq.to.split.revcomp = strsplit(revCompString(check.seq), split = "")[[1]]
          #handle problematic first and last cases
          shifts$seq[1] = paste(seq.to.split.plus, collapse = "")
          shifts$seq[(seq.len)] = paste(c(seq.to.split.plus[seq.len], seq.to.split.plus[1 : (seq.len - 1)]), collapse = "")
          shifts$seq[seq.len  + 1] = paste(seq.to.split.revcomp, collapse = "")
          shifts$seq[(seq.len * 2)] = paste(c(seq.to.split.revcomp[seq.len], seq.to.split.revcomp[1 : (seq.len - 1)]), collapse = "")
          #do forward strand
          for(ii in 2 : (seq.len - 1))
          {
            shifts$seq[ii] = paste(c(seq.to.split.plus[ii : seq.len], seq.to.split.plus[1 : (ii - 1)]), collapse = "")
          }
          #do reverse strand
          for(ii in 2 : (seq.len - 1))
          {
            shifts$seq[ii + seq.len] = paste(c(seq.to.split.revcomp[ii : seq.len], seq.to.split.revcomp[1 : (ii - 1)]), collapse = "")
          }
          #align sequence.template to the shifts
          for(ii in 1 : nrow(shifts))
          {
            if(ii%%100 == 0)
            {
              print(paste("Aligning shfit ", ii, " out of ", nrow(shifts), " from region ", i, paste = ""))
            }
            #TODO choose an alignment method that gives a precise score and align and save the score in the table
            shifts$score[ii] =  pairwiseAlignment(pattern = shifts$seq[ii], subject = sequence.template$seq[highest], type = "global", scoreOnly = TRUE)
          }
          #return the highest score shift as a primary.shifted[i]
          regions.data.frame$consensus.primary[i] = shifts$seq[which.max(shifts$score)]
          regions.data.frame$class[i] = sequence.template$name[highest]
        }
      }
    }
  }
  
  
  # 4. Refine each consensus by mapping again, aligning with mafft and extracting consensus 
  
  if(temp.folder == "")
  {
    stop("no temp folder")
  }
  setwd(temp.folder)
  
  if(!dir.exists(assemblyName))
  {
    dir.create(assemblyName)
  }
  if(!dir.exists(paste(assemblyName, "/frame4", sep ="")))
  {
    dir.create(paste(assemblyName, "/frame4", sep =""))
  }
  if(!dir.exists(paste(assemblyName, "/frame4/inputs", sep = "")))
  {
    dir.create(paste(assemblyName, "/frame4/inputs", sep = ""))
  }
  if(!dir.exists(paste(assemblyName, "/frame4/outputs", sep = "")))
  {
    dir.create(paste(assemblyName, "/frame4/outputs", sep = ""))
  }
  
  consensus.secondary = vector(mode = "character", length = nrow(regions.data.frame))
  repeats.identified = vector(mode = "numeric", length = nrow(regions.data.frame))
  
  
  for(i in 1 : nrow(regions.data.frame))
  {
    if(regions.data.frame$consensus.primary[i] != "none_identified")
    {
      print(paste("on ",regions.data.frame$fasta.name[i], ", generating consensus for region ", i, "/", nrow(regions.data.frame), " window size: ", regions.data.frame$end[i] - regions.data.frame$start[i], sep = ""))
      seqC = str_sub(DNA.sequence, regions.data.frame$start[i], regions.data.frame$end[i])
      N = regions.data.frame$most.freq.value.N[i]
      maxMis = N %/% 3
      if(maxMis > 100)
      {
        maxMis = 100
      }
      matchPlus = matchPattern(pattern = regions.data.frame$consensus.primary[i], subject = seqC, max.mismatch = maxMis, with.indels = TRUE)
      tempP = as.data.frame(matchPlus)
      tempP$start  = start(matchPlus)
      tempP$end = end(matchPlus)
      if(nrow(tempP) > 0)
      {
        tempP$strand = "+"
      }
      matchMinus = matchPattern(pattern = revCompString(regions.data.frame$consensus.primary[i]), subject = seqC, max.mismatch = maxMis, with.indels = TRUE)
      tempM = as.data.frame(matchMinus)
      tempM$start  = start(matchMinus)
      tempM$end = end(matchMinus)
      if(nrow(tempM) > 0)
      {
        tempM$strand = "-"
      }
      match = rbind(tempP, tempM)
      match = match[order(match$start),]
      
      ########## handle overlaps
      
      if(nrow(match) > 0)
      {
        iii = nrow(match)
        
        while(iii > 1)
        {
          overlap = match$start[iii] - match$end[iii - 1] - 1
          
          if(overlap > 0 & overlap <= 10)    #if ovarlap is positive, it means its a gap and small gap will be bridged by adding up to 10 nt to the end of the previous seq
          {
            match$end[iii-1] = match$end[iii-1] + overlap
            
          } else if(overlap < 0 & overlap >= -10) # if overlap is small, seq will be divided evenly between two repeats
          {
            overlap = -overlap
            b = overlap %/% 2
            c = overlap - b
            match$start[iii] = match$start[iii] + b
            match$end[iii - 1] = match$end[iii - 1] - c
            
          } else if(overlap < -10 & overlap > -500) # if overlap is big, the next seq will be removed completely
          {
            match = match[-c(iii),] 
          }
          iii = iii - 1
        }
      }
      
      if(nrow(match) > 0)
      {
        for(ii in 1 : nrow(match))
        {
          match$seq[ii] = substr(seqC, start = match$start[ii], stop = match$end[ii])
          {
            if(match$strand[ii] == "-")
            {
              match$seq[ii] = revCompString(match$seq[ii])
            }
          }
        }
        match$width = nchar(match$seq)
        match = match[match$end > match$start,]
        match = match[match$width > 0,]
        match$class = regions.data.frame$class[i]
        match$region.name = paste(assemblyName, fasta.name, sep = "_")
        
      }
      
      
      if(nrow(match) > 0)
      {
        write.csv(x = match, file = paste(assemblyName, "/frame4", "/Repeats_", regions.data.frame$name[i], "_", regions.data.frame$start[i], "_", regions.data.frame$end[i], ".csv", sep = ""))
        #export sequences in fasta
        #for(ii in 1 : nrow(match))
        #{
        #  write.fasta(sequences = match$seq[ii], names = paste(assemblyName, ".primary.extract.", regions.data.frame$index[i], ".", regions.data.frame$name[i], ".", ii, sep = ""), 
        #              file.out = paste(assemblyName, "/frame4/inputs/primary.extract", ".", regions.data.frame$index[i], ".", regions.data.frame$name[i], ".", "fasta", sep = ""), open = "a", as.string = TRUE)
        #}
        write.fasta(sequences = str_split(match$seq, pattern = ""), names = paste(assemblyName, ".primary.extract.", regions.data.frame$index[i], ".", regions.data.frame$name[i], ".", match$start, sep = ""), 
                    file.out = paste(assemblyName, "/frame4/inputs/primary.extract", ".", regions.data.frame$index[i], ".", regions.data.frame$name[i], ".", "fasta", sep = ""), open = "w")
        
        
        #mafft
        
        input = paste(assemblyName, "/frame4/inputs/primary.extract", ".", regions.data.frame$index[i], ".", regions.data.frame$name[i], ".", "fasta", sep = "")
        output = paste(assemblyName, "/frame4/outputs/primary.extract", ".", regions.data.frame$index[i], ".", regions.data.frame$name[i], ".", "aligned.fasta", sep = "")
        if(Sys.info()['sysname'] == "Linux")
        {
          system(paste("mafft --quiet --retree 2 --inputorder ", input, " > ", output, sep = ""), intern = TRUE)
        } else if(Sys.info()['sysname'] == "Windows")
        {
          shell(paste("mafft --quiet --retree 2 --inputorder ", input, " > ", output, sep = ""), intern = TRUE)
        }
        
        alignment = read.alignment(output, format = "FASTA", forceToLower = FALSE)
        consensus = consensus(alignment, threshold = 0.3)
        consensus = toupper(consensus[consensus != "-"])
        consensus = paste(consensus, collapse = "")
        repeats.identified[i] = length(alignment$seq)
        
        consensus.secondary[i] = consensus
        remove(consensus)
      }
    } else
    {
      consensus.secondary[i] = "none_identified"
      repeats.identified[i] = 0
    }
  }
  if(nrow(regions.data.frame) > 0)
  {
    regions.data.frame = cbind(regions.data.frame, consensus.secondary, repeats.identified)
  }
  gc()
  return(regions.data.frame)
}

draw.scaffold.repeat.plots = function(temp.folder = "", assemblyName = "", fastaDirectory = "", only.Chr1_5 = TRUE, single.pngs = FALSE)
{
  fasta = read.fasta(fastaDirectory, as.string = TRUE)
  if(temp.folder == "")
  {
    stop("no temp folder")
  }
  setwd(temp.folder)
  
  if(!dir.exists(assemblyName))
  {
    dir.create(assemblyName)
  }
  
  if(!dir.exists(paste(temp.folder, "/", assemblyName, "/plots", sep = "")))
  {
    dir.create(paste(temp.folder, "/", assemblyName, "/plots", sep = ""))
  }
  
  
  
  files = list.files(path = paste(temp.folder, "/", assemblyName, "/frame4", sep = ""), pattern = "Repeats*")
  
  if(only.Chr1_5 == TRUE)
  {
    a = as.data.frame(str_split(files, pattern = "_", simplify = TRUE))
    files = files[which(a$V2 == "Chr1" | a$V2 == "Chr2" | a$V2 == "Chr3" | a$V2 == "Chr4" | a$V2 == "Chr5")]
  }
  
  rbPal = colorRampPalette(c('red','blue'))
  
  if(single.pngs == FALSE)
  {
    #save as a single pdf
    i = 1
    pdf(file = paste(temp.folder, "/", assemblyName, "/plots/", assemblyName, ".pdf", sep = ""), width = 50, height = 10, pointsize = 30)
    while(i < length(files))
    {
      k = i
      scaffold = str_split(files[i], pattern = "_")[[1]]
      scaffold = paste(scaffold[2 : (length(scaffold) - 2)], collapse = "_")
      fragment.length = as.numeric(nchar(fasta[which(names(fasta) == scaffold)]))
      if(length(fragment.length) != 0)
      {
        plot(0,0, type="p", xlab="repeat start position", ylab="repeat size", xlim=c(0, fragment.length), ylim=c(0, 600), main = scaffold, cex = 0)
        
        while(k <= length(files) && paste((str_split(files[k], pattern = "_")[[1]])[2 : (length((str_split(files[k], pattern = "_")[[1]])) - 2)], collapse = "_") == scaffold)
        {
          repeats = read.csv(file = paste(temp.folder, "/", assemblyName, "/frame4/", files[k], sep = ""), header = TRUE)
          repeats$start = repeats$start + as.numeric(str_split(files[k], pattern = "_")[[1]][length(str_split(files[k], pattern = "_")[[1]]) - 1])
          repeats$col = rbPal(600)[repeats$width]
          points(repeats$start,repeats$width, col = repeats$col, pch = ".", cex = 7)
          remove(repeats)
          k = k + 1
        }
      }
      
      i = k
    }
    dev.off()
  }
  if(single.pngs == TRUE)
  {
    i = 1
    while(i < length(files))
    {
      k = i
      scaffold = str_split(files[i], pattern = "_")[[1]]
      scaffold = paste(scaffold[2 : (length(scaffold) - 2)], collapse = "_")
      fragment.length = as.numeric(nchar(fasta[which(names(fasta) == scaffold)]))
      
      if(length(fragment.length) > 0)
      {
        png(filename = paste(temp.folder, "/", assemblyName, "/plots/", assemblyName, "_", i, ".png", sep = ""), width = 4000, height = 1000, pointsize = 45)
        plot(0,0, type="p", xlab="repeat start position", ylab="repeat size", xlim=c(0, fragment.length), ylim=c(0, 600), main = scaffold, cex = 0)
        
        while(k <= length(files) && paste((str_split(files[k], pattern = "_")[[1]])[2 : (length((str_split(files[k], pattern = "_")[[1]])) - 2)], collapse = "_") == scaffold)
        {
          repeats = read.csv(file = paste(temp.folder, "/", assemblyName, "/frame4/", files[k], sep = ""), header = TRUE)
          repeats$start = repeats$start + as.numeric(str_split(files[k], pattern = "_")[[1]][length(str_split(files[k], pattern = "_")[[1]]) - 1])
          repeats$col = rbPal(600)[repeats$width]
          points(repeats$start,repeats$width, col = repeats$col, pch = ".", cex = 7)
          remove(repeats)
          k = k + 1
        }
        dev.off()
      }
      i = k
    }
  }
}

extract.all.repeats = function(temp.folder = "", assemblyName = "")
{
  if(temp.folder == "")
  {
    stop("no temp folder")
  }
  print(paste("start extract repeats from ", assemblyName, sep = ""))
  setwd(temp.folder)
  start = vector(mode = "numeric", length = 0)
  end = vector(mode = "numeric", length = 0)
  width = vector(mode = "numeric", length = 0)
  seq = vector(mode = "character", length = 0)
  region.name = vector(mode = "character", length = 0)
  strand = vector(mode = "character", length = 0)
  class =  vector(mode = "numeric", length = 0)
  
  repeats.all = data.frame(start, end, width, seq, strand, class, region.name)
  remove(start, end, strand, seq, width, class, region.name)
  
  files = list.files(path = paste(assemblyName, "/frame4", sep = ""), pattern = "Repeats_*")
  if(length(files) == 0)
  {
    return(0)
  }
  for(i in 1 : length(files))
  {
    region = read.csv(file = paste(assemblyName, "/frame4/", files[i], sep = ""), header = TRUE)
    if(Sys.info()['sysname'] == "Linux")
    {
      region = region[,-c(1)]
    } else if(Sys.info()['sysname'] == "Windows")
    {
      region = region[,-c(1,2)]
    }
    
    
    #region = region[,c(1,3,4,7,6,5)]
    region.coor = strsplit(files[i], split = "_")[[1]]
    region.start = as.numeric(region.coor[length(region.coor)-1])
    region.end = as.numeric(strsplit(region.coor[length(region.coor)], split = "\\.")[[1]][1])
    region$start = region$start + region.start
    region$end = region$end + region.start
    
    repeats.all = rbind(repeats.all, region)
  }
  repeats.all = repeats.all[order(repeats.all$start),]
  repeats.all = repeats.all[order(repeats.all$region.name),]
  write.csv(x = repeats.all, file = paste(assemblyName, "/all.repeats.from.", assemblyName, ".csv", sep = ""), row.names = FALSE)
}




#####

#####
#RUN#
#####

if(test.sequence == TRUE)
{
  if(Sys.info()['sysname'] == "Linux")
  {
    outputs.directory = "~/Arabidopsis/output" 
    genomes.directory = "~/Arabidopsis/genomes" 
    #sequence.templates = "~/Arabidopsis/sequence.template.csv"
  } else if(Sys.info()['sysname'] == "Windows")
  {
    outputs.directory = "C:/Users/wlodz/Desktop/RepeatIdentifyTests/output" 
    genomes.directory = "C:/Users/wlodz/Desktop/RepeatIdentifyTests/genomes" 
    #sequence.templates = "C:/Users/wlodz/Desktop/RepeatIdentifyTests/sequence.template.csv"
  }
}
if(!is.na(sequence.templates))
{
  sequence.templates = read.csv(file = sequence.templates)
  sequence.templates = sequence.templates[sequence.templates$group == "ath",]
}


if(Sys.info()['sysname'] == "Linux")
{
  genome.names = system(paste("cd ", genomes.directory, " && ls -d *.fa*", sep = ""), intern = TRUE)
} else if(Sys.info()['sysname'] == "Windows")
{
  genome.names = shell(paste("cd ", genomes.directory, " && dir /b | findstr .fa", sep = ""), intern = TRUE)
}

sequences = data.frame(index = numeric(), assembly.name = character(), fasta.name = character(), sequence = character())
index = 1
for(i in 1 : length(genome.names))
{
  print(paste("reading file ", genome.names[i], sep = ""))
  fasta = read.fasta(file = paste(genomes.directory, genome.names[i], sep = "/"), as.string = TRUE)
  for(j in 1 : length(fasta))
  {
    print(paste("reading sequence in ", genome.names[i], " ", j, " out of ", length(fasta), sep = ""))
    sequences = rbind(sequences, data.frame(index = index, assembly.name = genome.names[i], fasta.name = names(fasta)[j], sequence = toupper(fasta[[j]][1])))
    index = index + 1
  }
}
if(Sys.info()['sysname'] == "Linux")
{
  if(set.no.of.cores != 0)
  {
    registerDoParallel(cores = set.no.of.cores)
  } else {
    registerDoParallel(cores = nrow(sequences))
  }
} else if(Sys.info()['sysname'] == "Windows")
{
}

#####
if(Sys.info()['sysname'] == "Linux")
{
  print("Currently registered parallel backend name, version and cores")
  print(getDoParName())
  print(getDoParVersion())
  print(getDoParWorkers())
  foreach(i = 1 : nrow(sequences)) %dopar% {
    print(paste("Assembly ", sequences$assembly.name[i], " chromosome ", sequences$fasta.name[i], " started", sep = ""))
    identified = RepeatIdentifier(DNA.sequence = sequences$sequence[i], assemblyName = sequences$assembly.name[i], fasta.name = sequences$fasta.name[i], 
                                  kmer = set.kmer, window = 1000, threshold = set.threshold, mask.small.regions = 1500, mask.small.repeats = 4,
                                  max.repeat.size = set.max.repeat.size,
                                  tests = 3, temp.folder = outputs.directory, sequence.template = sequence.templates)
    if((identified != 0)[1])
    {
      print(paste("Assembly ", sequences$assembly.name[i], " chromosome ", sequences$fasta.name[i], " finished", sep = ""))
      identified$chr.no = i
      write.csv(x = identified, file = paste(outputs.directory, "/", sequences$assembly.name[i], "/Chromosome.", sequences$fasta.name[i], ".csv", sep = ""), quote = FALSE)
    }
  } 
} else if(Sys.info()['sysname'] == "Windows")
{
  for(i in 1 : nrow(sequences)){
    print(paste("Assembly ", sequences$assembly.name[i], " chromosome ", sequences$fasta.name[i], " started", sep = ""))
    identified = RepeatIdentifier(DNA.sequence = sequences$sequence[i], assemblyName = sequences$assembly.name[i], fasta.name = sequences$fasta.name[i], 
                                  kmer = set.kmer, window = 1000, threshold = set.threshold, mask.small.regions = 1500, mask.small.repeats = 4,
                                  max.repeat.size = set.max.repeat.size,
                                  tests = 3, temp.folder = outputs.directory, sequence.template = sequence.templates)
    if((identified != 0)[1])
    {
      print(paste("Assembly ", sequences$assembly.name[i], " chromosome ", sequences$fasta.name[i], " finished", sep = ""))
      identified$chr.no = i
      write.csv(x = identified, file = paste(outputs.directory, "/", sequences$assembly.name[i], "/Chromosome.", sequences$fasta.name[i], ".csv", sep = ""), quote = FALSE)
    }
  } 
}
assemblies = unique(sequences$assembly.name)

for(i in 1 : length(assemblies))
{
  print(i)
  if(Sys.info()['sysname'] == "Linux")
  {
    repeat.files = system(paste("cd ", outputs.directory, "/", genome.names[i], " && ls -d Chromosome.*", sep = ""), intern = TRUE)
    print(repeat.files)
    regions = NULL
    for(j in 1 : length(repeat.files))
    {
      regions = rbind(regions, read.csv(file = paste(outputs.directory, "/", genome.names[i], "/", repeat.files[j], sep = "")))
    }
    regions = regions[,-c(1)]
    print("writing regions data frame")
    write.csv(x = regions, file = paste(outputs.directory, "/", genome.names[i], "/Summary.of.repetitive.regions.", genome.names[i], ".csv", sep = ""), quote = FALSE)
    
    print("drawing plots")
    draw.scaffold.repeat.plots(temp.folder = outputs.directory, assemblyName = assemblies[i], 
                               fastaDirectory = paste(genomes.directory, "/", assemblies[i], sep = ""), only.Chr1_5 = FALSE, single.pngs = TRUE)
    
    extract.all.repeats(temp.folder = outputs.directory, assemblyName = assemblies[i])
    print("all done for this genome")
  } else if(Sys.info()['sysname'] == "Windows")
  {
    repeat.files = shell(paste("cd ", outputs.directory, "/", genome.names[i], " && dir /b | findstr Chromosome.", sep = ""), intern = TRUE)
    
    regions = NULL
    for(j in 1 : length(repeat.files))
    {
      regions = rbind(regions, read.csv(file = paste(outputs.directory, "/", genome.names[i], "/", repeat.files[j], sep = "")))
    }
    regions = regions[,-c(1,2)]
    print("writing regions data frame")
    write.csv(x = regions, file = paste(outputs.directory, "/", genome.names[i], "/Summary.of.repetitive.regions.", genome.names[i], ".csv", sep = ""), quote = FALSE)
    
    print("drawing plots")
    draw.scaffold.repeat.plots(temp.folder = outputs.directory, assemblyName = assemblies[i], 
                               fastaDirectory = paste(genomes.directory, "/", assemblies[i], sep = ""), only.Chr1_5 = FALSE, single.pngs = TRUE)
    
    extract.all.repeats(temp.folder = outputs.directory, assemblyName = assemblies[i])
  }
}






