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

#Parallelisation#
  #it's based on amount of .fasta files provided as an input, to use more CPUs when running a single genome, it can be split into chromosomes
#if left at 0, R will try to use the same amount or cores as the amount of .fasta files in the input folder
set.no.of.cores = 0 

#Run settings#
set.kmer = 12 # kmer size used for initial identification of repetitive regions
set.threshold = 10 # window repetitiveness score (0-100) threshold
set.max.repeat.size = 500 # max size of repeats to be identified
outputs.directory = "" # output folder. example "~/Arabidopsis/output"
genomes.directory = "" # folder with .fasta inputs. example "~/Arabidopsis/genomes"
sequence.templates = NA # used to match identified repeats to templates so they are globally in the same frame (same start position), set as NA to skip
  ##example: sequence.templates = read.csv(file = "~/sequence.template.csv")
  ##format: 
  ##seq,name,length,group
  ##AGTATA,CEN180,180,ath

#MAFFT#
#without mafft (https://mafft.cbrc.jp/alignment/software/) installed in the environment only regions containing repeats will be identified
#tested with mafft-7.475-win64-signed
skip.mafft = FALSE



test.sequence = TRUE # keep FALSE, otherwise outputs and genomes directory paths will be overwritten
###########
#Functions#
###########

revCompString = function(DNAstr) {return(toupper(toString(reverseComplement(DNAString(DNAstr))))) }

# 1. Identify repetitive regions
identify.repetitive.regions = function(fastaDirectory = "", kmer = 12, window = 1000, threshold = 10, filter.small = 0)
{
  time = Sys.time()
  
  fasta = read.fasta(fastaDirectory, as.string = TRUE)
  
  repetitive.windows = data.frame(index = integer(), fasta.name = character(), start = integer(), end = integer(), average.score = double(), window.sequence = character(), most.freq.value.N = integer())
  
  index = 0
  
  for(i in 1 : length(fasta))   # for every fasta sequence in the file
  {
    
    seqA = fasta[[i]]
    seqA = toupper(seqA)
    seq.length = nchar(seqA)
    scores = vector(mode = "numeric", length = seq.length %/% window + 1)
    
    for(ii in 1 : (length(scores) - 1))  # for every window of size window but not the last one
    {
      if(ii%%window == 0)
      {
        print(paste("sequence ", i, "/", length(fasta), ", window ", ii, "/", length(scores), sep = ""))
      }
      window.start = (ii - 1) * window
      check.window = str_sub(seqA, window.start, (window.start + window - 1))
      for(iii in 1 : (window - kmer))
      {
        check.string = str_sub(string = check.window, start = iii, end = (iii + kmer - 1))
        count = 1 * ((str_count(string = check.string, pattern = "N") == 0) & str_count(string = check.window, pattern = check.string) > 1)
        scores[ii] = scores[ii] + count
      }
    }
    ii = length(scores) # for the last window
    print(paste("sequence ", i, "/", length(fasta), ", window ", ii, "/", length(scores), sep = ""))
    window.start = (ii - 1) * window
    check.window = str_sub(seqA, window.start)
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
        name = names(fasta)[i]
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
        window.sequence = str_sub(seqA, start, end)
        most.freq.value.N = 0
        temp = data.frame(index, name, start, end, ave.score, window.sequence, most.freq.value.N)
        repetitive.windows = rbind(repetitive.windows, temp)
      }
      ii = ii + 1
    }
    
    repetitive.windows = repetitive.windows[(repetitive.windows$end - repetitive.windows$start) > filter.small,]
    
    print(paste("Step 1 (identify windows) done in ", (Sys.time() - time), sep = ""))
  }
  gc()
  return(repetitive.windows)
}

# 2. Calculate closest identical k-mer distances
closest.identical.kmer.distances = function(regions.data.frame, kmer = 12, reach = 500, filter.smaller.N = 2, plot = FALSE)
{
  for(i in 1 : nrow(regions.data.frame))
  {
    seqB = regions.data.frame$window.sequence[i]
    print(paste("calculating distances on region ", i, "/", nrow(regions.data.frame), " region size: ", nchar(seqB), sep = ""))
    distance = vector(mode = "numeric", length =  nchar(seqB))
    for(ii in 1 : (nchar(seqB) - kmer + 1))
    {
      kmer.pattern = str_sub(seqB, ii, ii + kmer - 1)
      window.string = str_sub(seqB, ii + 1, ii + reach)
      distance[ii] = str_locate(string = window.string, pattern = kmer.pattern)[[1]]
    }
    distance = distance[!is.na(distance)]
    distance = distance[distance > 0]
    if(length(distance) > 0)
    {
      a = hist(distance, breaks = max(distance), plot = FALSE)
      regions.data.frame$most.freq.value.N[i] = a$breaks[which.max(a$counts) + 1]
      if(plot == TRUE)
      {
        png(filename = paste("", regions.data.frame$name[i], ".", regions.data.frame$index[i], ".N=", regions.data.frame$most.freq.value.N[i], ".png", sep = ""), width = 750, height = 500, pointsize = 25)
        hist(distance, breaks = max(distance), main = paste("Histogram of", regions.data.frame$name[i], "index", regions.data.frame$index[i], "N =", regions.data.frame$most.freq.value.N[i], sep = " "))
        dev.off()
      }
    }
  }
  ##regions.data.frame = regions.data.frame[regions.data.frame$most.freq.value.N >= filter.smaller.N,]
  gc()
  return(regions.data.frame)
}

# 3. Test random N-size mers, handle overlaps, choose the best one, extract repeats for mafft alignment
test.random.Nmers = function(regions.data.frame, tests = 5, assemblyName, temp.folder = "", remove.temp = TRUE)
{
  #regions.data.frame$full.name[i]
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
    print(paste("testing region ", i, "/", nrow(regions.data.frame), " window size: ", regions.data.frame$end[i] - regions.data.frame$start[i], sep = ""))
    seqC = regions.data.frame$window.sequence[i]
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
        if(temp$start[iii] - temp$end[iii - 1] < 0)
        {
          temp = temp[-iii,]
        }
        iii = iii - 1
      }
      write.csv(temp, file = paste(temp.folder, "/", assemblyName, "/frame2/", regions.data.frame$full.name[i], "_", i, "_", ii, ".csv", sep = ""), row.names = FALSE)
      random.sequence.scores[ii] = sum(temp$end - temp$start)
    }
    #mafft
    winner = read.csv(paste(temp.folder, "/", assemblyName, "/frame2/", regions.data.frame$full.name[i], "_",  i, "_", which.max(random.sequence.scores), ".csv", sep = ""))
    for(ii in 1 : nrow(winner))
    {
      if(ncol(winner) == 5)
      {
        write.fasta(file.out = paste(assemblyName, "/", i, "_", "pre.extract", regions.data.frame$index[i], ".", regions.data.frame$name[i], ".fasta", sep = ""), 
                    names = paste(winner$start[ii], winner$end[ii], sep = "_"), sequences = winner[ii,4], open = "a", as.string = TRUE)
      } else if(ncol(winner) == 4){
        write.fasta(file.out = paste(assemblyName, "/", i, "_", "pre.extract", regions.data.frame$index[i], ".", regions.data.frame$name[i], ".fasta", sep = ""), 
                    names = paste(winner$start[ii], winner$end[ii], sep = "_"), sequences = winner$x[ii], open = "a", as.string = TRUE)
      }
    }
    
    input = paste(assemblyName, "/", i, "_", "pre.extract", regions.data.frame$index[i], ".", regions.data.frame$name[i], ".fasta", sep = "")
    output = paste(assemblyName, "/", i, "_", "pre.extract", regions.data.frame$index[i], ".", regions.data.frame$name[i], ".aligned.fasta", sep = "")
    
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
  }
  if(remove.temp)
  {
    unlink(assemblyName, recursive = TRUE)
  }
  print(paste("finished test random Nmers of", regions.data.frame$full.name[i], sep = " "))
  gc()
  return(cbind(regions.data.frame, consensus.primary, consensus.count))
}


#kamers

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

# 3.5 Shift the primary consensus sequences to match a given sequence, only those regions with N within a given range
shift.primary.consensus = function(regions.data.frame, sequence.template)
{
  for(i in 1 : nrow(regions.data.frame))
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
    }
  }
  return(cbind(regions.data.frame))
}

# 4. Refine each consensus by mapping again, aligning with mafft and extracting consensus 
#TODO: check for mafft version
#TODO: change hwere input and output mafft files are created to a different temp directory
#TODO: mafft.directory not needed
#TODO: check and remove .bat from, mafft
generate.secondary.consensus = function(regions.data.frame, assemblyName = "", temp.folder = "", temp.remove = FALSE)
{
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
    
    print(paste("generating consensus for region ", i, "/", nrow(regions.data.frame), " window size: ", regions.data.frame$end[i] - regions.data.frame$start[i], sep = ""))
    seqC = regions.data.frame$window.sequence[i]
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
    matchMinus = matchPattern(pattern = toupper(revCompString(regions.data.frame$consensus.primary[i])), subject = seqC, max.mismatch = maxMis, with.indels = TRUE)
    tempM = as.data.frame(matchMinus)
    tempM$start  = start(matchMinus)
    tempM$end = end(matchMinus)
    if(nrow(tempM) > 0)
    {
      tempM$strand = "-"
    }
    match = rbind(tempP, tempM)
    match = match[order(match$start),]
    
    ##########
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
    
    for(ii in 1 : nrow(match))
    {
      match$seq[ii] = substr(regions.data.frame$window.sequence[i], start = match$start[ii], stop = match$end[ii])
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
    match$region.name = regions.data.frame$full.name[i]
    
    
    
    if(nrow(match) > 0)
    {
      write.csv(x = match, file = paste(assemblyName, "/frame4", "/Repeats_", regions.data.frame$name[i], "_", regions.data.frame$start[i], "_", regions.data.frame$end[i], ".csv", sep = ""))
      #export sequences in fasta
      for(ii in 1 : nrow(match))
      {
        write.fasta(sequences = match$seq[ii], names = paste(assemblyName, ".primary.extract.", regions.data.frame$index[i], ".", regions.data.frame$name[i], ".", ii, sep = ""), 
                    file.out = paste(assemblyName, "/frame4/inputs/primary.extract", ".", regions.data.frame$index[i], ".", regions.data.frame$name[i], ".", "fasta", sep = ""), open = "a", as.string = TRUE)
      }
      
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
    if(temp.remove)
    {
      unlink(paste(assemblyName, "/frame4/inputs", sep = ""), recursive = TRUE)
      unlink(paste(assemblyName, "/frame4/outputs", sep = ""), recursive = TRUE)
    }
  }
  gc()
  return(cbind(regions.data.frame, consensus.secondary, repeats.identified))
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
    print(paste("Reading file ", i, " out of ", length(files), paste = ""))
    region = read.csv(file = paste(assemblyName, "/frame4/", files[i], sep = ""), header = TRUE)
    region = region[,-c(1)]
    
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


identify.and.extract.repeats = function(temp.folder = "", genome.dir = "", genome = "", sequence.templates, kmer = 10, window = 1000, threshold = 10, filter.small = 1000, reach = 200, 
                                        filter.smaller.N = 3, plot.regions.scores = FALSE, random.Nmer.tests = 3, only.Chr1_5 = FALSE, single.pngs = TRUE, debug.text = FALSE)
{
  print(paste("Start file ", genome, sep = ""))
  fasta = paste(genome.dir, genome, sep ="/")
  assemblyName = strsplit(genome, split = "[.]")[[1]]
  assemblyName = str_c(assemblyName[1:(length(assemblyName) - 1)], collapse = ".")
  a = Sys.time()
  frame1 = identify.repetitive.regions(fastaDirectory = fasta, kmer = kmer, window = window, threshold = threshold, filter.small = filter.small)
  frame1$full.name = assemblyName
  if(debug.text) {print("frame1 done")}
  frame2 = closest.identical.kmer.distances(regions.data.frame = frame1, kmer = kmer, reach = reach, filter.smaller.N = filter.smaller.N, plot = plot.regions.scores)
  if(debug.text) {print("frame2 done")}
  if(skip.mafft) 
  {
    frame7 = frame2[,-c(6)]
    write.csv(x = frame7, file = paste(temp.folder, "/", assemblyName, "/Summary.of.repetitive.regions.", assemblyName, ".csv", sep = ""), quote = FALSE)
    remove(frame1,frame2,frame7,fasta)
    b = difftime(Sys.time(), a, units = "hours")
    print(paste("Genome ", assemblyName, " finished in ", round(b, digits = 3), " hours", sep = ""))
    remove(assemblyName,a,b)
  } else {
    frame3 = test.random.Nmers(regions.data.frame = frame2, tests = random.Nmer.tests, assemblyName = assemblyName, temp.folder = temp.folder)
    if(!is.na(sequence.templates))
    {
      frame3.5 = shift.primary.consensus(regions.data.frame = frame3, sequence.template = sequence.templates)
      frame4 = generate.secondary.consensus(regions.data.frame = frame3.5, assemblyName = assemblyName, temp.folder = temp.folder)
      remove(frame3)
    } else{
      frame4 = generate.secondary.consensus(regions.data.frame = frame3, assemblyName = assemblyName, temp.folder = temp.folder)
    }
    
    frame7 = frame4[,-c(6)]
    write.csv(x = frame7, file = paste(temp.folder, "/", assemblyName, "/Summary.of.repetitive.regions.", assemblyName, ".csv", sep = ""), quote = FALSE)
    
    draw.scaffold.repeat.plots(temp.folder = temp.folder, assemblyName = assemblyName,fastaDirectory = fasta,only.Chr1_5 = only.Chr1_5,single.pngs = single.pngs)
    extract.all.repeats(temp.folder = temp.folder, assemblyName = assemblyName)
    if(debug.text) {print("extracted repeats done")}
    b = difftime(Sys.time(), a, units = "hours")
    remove(frame1,frame2,frame3,frame4,frame7,fasta)
    print(paste("Genome ", assemblyName, " finished in ", round(b, digits = 3), " hours", sep = ""))
    remove(assemblyName,a,b)
  }
}


#####

#####
#RUN#
#####

if(test.sequence == TRUE)
{
  if(Sys.info()['sysname'] == "Linux")
  {
    outputs.directory = "~/Arabidopsis/output" # output folder. example "~/Arabidopsis/output"
    genomes.directory = "~/Arabidopsis/genomes" # folder with .fasta inputs. example "~/Arabidopsis/genomes"
  } else if(Sys.info()['sysname'] == "Windows")
  {
    outputs.directory = "C:/Users/wlodz/Desktop/RepeatIdentifyTests/output" # output folder. example "~/Arabidopsis/output"
    genomes.directory = "C:/Users/wlodz/Desktop/RepeatIdentifyTests/genomes" # folder with .fasta inputs. example "~/Arabidopsis/genomes"
  }
}

if(Sys.info()['sysname'] == "Linux")
{
  genome.names = system(paste("cd ", genomes.directory, " && ls -d *.fa*", sep = ""), intern = TRUE)
} else if(Sys.info()['sysname'] == "Windows")
{
  genome.names = shell(paste("cd ", genomes.directory, " && dir /b | findstr .fa", sep = ""), intern = TRUE)
}
if(set.no.of.cores != 0)
{
  registerDoParallel(cores = set.no.of.cores)
} else {
  registerDoParallel(cores = length(genome.names))
}

if(Sys.info()['sysname'] == "Linux")
{
  print("Currently registered parallel backend name, version and cores")
  print(getDoParName())
  print(getDoParVersion())
  print(getDoParWorkers())
  print("lets go")
  foreach(i = 1 : length(genome.names)) %dopar% {
    identify.and.extract.repeats(temp.folder = outputs.directory, genome.dir = genomes.directory, genome = genome.names[i], 
                                 sequence.templates = sequence.templates, kmer = set.kmer, threshold = set.threshold, reach = set.max.repeat.size, debug.text = FALSE)
  } 
} else if(Sys.info()['sysname'] == "Windows")
{
  
  genome.names = shell(paste("cd ", genomes.directory, " && dir /b | findstr .fa", sep = ""), intern = TRUE)
  for(i in 1 : length(genome.names)){
    identify.and.extract.repeats(temp.folder = outputs.directory, genome.dir = genomes.directory, genome = genome.names[i], 
                                 sequence.templates = sequence.templates, kmer = set.kmer, threshold = set.threshold, reach = set.max.repeat.size, debug.text = FALSE)
  } 
}
#####

