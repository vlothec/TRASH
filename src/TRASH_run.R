#!/usr/bin/env Rscript
print("start")
thisFile = function() 
{
  cmd.Args = commandArgs(trailingOnly = FALSE)
  find.file = "--file="
  match = grep(find.file, cmd.Args)
  if (length(match) > 0) {
    return(normalizePath(sub(find.file, "", cmd.Args[match])))
  } else {
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}



######### 
installation.path = thisFile()
installation.path = strsplit(installation.path, split = "/")[[1]]
installation.path = paste(installation.path[1:(length(installation.path) - 2)], collapse = "/")
execution.path = getwd()
#attach R libraries


src.files = list.files(path = paste(installation.path, "/src", sep = ""), pattern = "fn_", full.names = TRUE)

print(Sys.time())
for(i in 1 : length(src.files))
{
  print(paste("attaching function from: ", src.files[i], sep = ""))
  source(src.files[i])
}
arguments = commandArgs(trailingOnly = TRUE)


#run settings
{
  PLOTTING.ONLY = FALSE
  sequence.templates = NA
  skip.repetitive.regions = FALSE
  change.lib.paths = TRUE
  cores.no = 1
  max.divergence.value = 5
  min.hor.value = 3
  skip.short.fasta.sequences = 0
  set.kmer = 10 # kmer size used for initial identification of repetitive regions
  set.threshold = 10 # window repetitiveness score (0-100) threshold
  set.max.repeat.size = 1000 # max size of repeats to be identified
  filter.small.regions = 2000 # repetitive windows smaller than this size will be removed (helps getting rid of regions with short duplications)
  filter.small.repeats = 4 # repetitive windows where dominant kmer distance is lower than this value will be removed (for example AT dinucleotide repeats)
  window.size = 1500 # how far apart kmers can be in the initial search for exact matches. No repeats larger than this will be identified
  MAX.CHROMOSOMES.TODO = 100000 
  hor.only = FALSE
  class.name.for.HOR = ""
  delete_temp_output = FALSE
  LIMIT.REPEATS.TO.ALIGN = 78000 #in base pairs
  simpleplot = FALSE
  random.seed = NULL
}
fasta.list = NULL
{
  if(length(arguments) == 0)
  {
    stop("At least one argument must be supplied", call. = FALSE)
  } else
  {
    for(i in 1 : length(arguments))
    {
      if(arguments[i] == "--def")
      {
        change.lib.paths = FALSE
      } else if(arguments[i] == "--skipr")
      {
        skip.repetitive.regions = TRUE
      } else if(arguments[i] == "--randomseed")
      {
        random.seed = as.numeric(arguments[i + 1])
      } else if(arguments[i] == "--horclass")
      {
        class.name.for.HOR = arguments[i + 1]
      } else if(arguments[i] == "--limrepno")
      {
        LIMIT.REPEATS.TO.ALIGN = as.numeric(arguments[i + 1])
      } else if(arguments[i] == "--horonly")
      {
        hor.only = TRUE
      } else if(arguments[i] == "--minhor")
      {
        min.hor.value = as.numeric(arguments[i + 1])
      } else if(arguments[i] == "--maxdiv")
      {
        max.divergence.value = as.numeric(arguments[i + 1])
      } else if(arguments[i] == "--maxchr")
      {
        MAX.CHROMOSOMES.TODO = as.numeric(arguments[i + 1])
      } else if(arguments[i] == "--k")
      {
        set.kmer = as.numeric(arguments[i + 1])
      } else if(arguments[i] == "--rmtemp")
      {
        delete_temp_output = TRUE
      } else if(arguments[i] == "--t")
      {
        set.threshold = as.numeric(arguments[i + 1])
      } else if(arguments[i] == "--m")
      {
        set.max.repeat.size = as.numeric(arguments[i + 1])
      } else if(arguments[i] == "--freg")
      {
        filter.small.regions = as.numeric(arguments[i + 1])
      } else if(arguments[i] == "--frep")
      {
        filter.small.repeats = as.numeric(arguments[i + 1])
      } else if(arguments[i] == "--o")
      {
        execution.path = arguments[i + 1]
      } else if(arguments[i] == "--win")
      {
        window.size = as.numeric(arguments[i + 1])
      } else if(arguments[i] == "--plotonly")
      {
        PLOTTING.ONLY = TRUE
      }
      else if(arguments[i] == "--seqt")
      {
        sequence.templates = arguments[i + 1]
      }
      else if(arguments[i] == "--simpleplot")
      {
        simpleplot = TRUE
      }
      else if(arguments[i] == "--par")
      {
        cores.no = as.numeric(arguments[i + 1])
      }
      if(grepl(".fa", arguments[i]) | grepl(".fna", arguments[i]))
      {
        if(file.exists(paste(execution.path, arguments[i], sep = "/")))
        {
          fasta.list = c(paste(execution.path, arguments[i], sep = "/"))
        } else if(file.exists(arguments[i]))
        {
          fasta.list = c(fasta.list, arguments[i])
        } else
        {
          stop(paste("Cannot find the", arguments[i], "fasta file specified", sep = ""), call. = FALSE)
        }
      }
    }
  }
  if(Sys.info()['sysname'] != "Linux")
  {
    print("Only Linux OS tested")
  } 
}

if(change.lib.paths)
{
  lib.path = paste(installation.path, "/libs", sep = "")
  .libPaths(lib.path)
}


print(.libPaths()[1])


suppressPackageStartupMessages(library(stringr, quietly = TRUE))
suppressPackageStartupMessages(library(stringdist, quietly = TRUE))
suppressPackageStartupMessages(library(base, quietly = TRUE))
suppressPackageStartupMessages(library(Biostrings, quietly = TRUE))
suppressPackageStartupMessages(library(seqinr, quietly = TRUE))
suppressPackageStartupMessages(library(doParallel, quietly = TRUE))
suppressPackageStartupMessages(library(circlize, quietly = TRUE))   #TODO add to the installer v0.4.15


if(PLOTTING.ONLY & hor.only)
{
  stop("Either only plot or only calculate HORs", call. = FALSE)
}

hor.c.script.path = paste(installation.path, "/src/HOR.V3.3", sep = "")


# read in fasta and setup parallel to the length of fasta file
print(fasta.list)
{
  fasta.sequence = NULL
  sequences = data.frame(file.name = "", fasta.name = "")
  for(i in 1 : length(fasta.list))
  {
    fasta = read.fasta(file = fasta.list[i], seqtype = "DNA", as.string = TRUE)
    if(length(fasta) < MAX.CHROMOSOMES.TODO)
    {
      skip.short.fasta.sequences = 0
    } else
    {
      if(length(nchar(fasta)[order(nchar(fasta), decreasing = TRUE)]) > MAX.CHROMOSOMES.TODO)
      {
        skip.short.fasta.sequences = nchar(fasta)[order(nchar(fasta), decreasing = TRUE)][MAX.CHROMOSOMES.TODO] 
      }
    }
    
    for(j in 1 : length(fasta))
    {
      if(nchar(fasta[j]) >= skip.short.fasta.sequences)
      {
        fasta.sequence = c(fasta.sequence,toupper(fasta[j]))
        sequences = rbind(sequences, data.frame(file.name = strsplit(fasta.list[i], split = "/")[[1]][length(strsplit(fasta.list[i], split = "/")[[1]])], fasta.name = names(fasta[j])))
        strsplit(fasta.list[i], split = "/")[[1]]
      }
    }
  }
  sequences = sequences[-1,]
  if(cores.no != 0)
  {
    registerDoParallel(cores = cores.no)
  } else
  {
    registerDoParallel(cores = nrow(sequences))
  }
  print("Currently registered parallel backend name, version and cores")
  print(getDoParName())
  print(getDoParVersion())
  print(getDoParWorkers())
}

if(!is.na(sequence.templates)) 
{
  if(file.exists(sequence.templates))
  {
    sequence.templates = read.csv(file = sequence.templates)
  } else
  {
    stop("Provided sequence template file does not exist")
  }
  names.check = c("name", "seq", "length") %in% names(sequence.templates)
  for(i in 1 : length(names.check))
  {
    if(!names.check[i])
    {
      print("Warning: sequence templates file might be incomplete")
    }
  }
}


if(hor.only)
{
  if(class.name.for.HOR == "")
  {
    stop("Repeat class not specified for HOR calculations", call. = FALSE)
  }
  
  print("calculating HORs")
  
  if(FALSE) #TODO bring back
  {
    #do HORs per chromosome
    foreach(i = 1 : length(fasta.sequence)) %dopar% {  
      #for(i in 1 : length(fasta.sequence)) {
      HOR.wrapper(threshold = 5, 
                  cutoff = 2, 
                  temp.folder = execution.path, 
                  assemblyName = sequences$file.name[i], 
                  chr.name = sequences$fasta.name[i],
                  mafft.bat.file = paste(installation.path, "/src/mafft-linux64/mafft.bat", sep = ""),
                  hor.c.script.path = hor.c.script.path,
                  class.name = class.name.for.HOR)
      gc()
    } 
  }
  
  
  for(i in 1 : length(fasta.sequence))
  {
    #do repetitiveness per chromosome
    calc.repetitiveness(temp.folder = execution.path, 
                        assemblyName = sequences$file.name[i], 
                        chr.name = sequences$fasta.name[i],
                        hor.c.script.path = hor.c.script.path,
                        class.name = class.name.for.HOR)
    #plot repetitiveness per chromosome
    plot.repetitiveness(temp.folder = execution.path, 
                                 assemblyName = sequences$file.name[i], 
                                 chr.name = sequences$fasta.name[i],
                                 hor.c.script.path = hor.c.script.path,
                                 class.name = class.name.for.HOR)
  }
  
  
  for(i in 1 : length(fasta.sequence))#### TODO remove this from here
  {
    #plot edit per chromosome
    plot.edit(temp.folder = execution.path, 
              assemblyName = sequences$file.name[i], 
              chr.name = sequences$fasta.name[i],
              hor.c.script.path = hor.c.script.path,
              class.name = class.name.for.HOR)
  }
  
  
  
  
  print("HORs finished)")
  quit()
  print("you can't see me")
}


if(PLOTTING.ONLY)
{
  #load
  for(i in 1 : length(fasta.list))
  {
    print(paste("plotting", sep = ""))
    
    outputs.directory <- paste(execution.path, sep = "")
    plot.min <- filter.small.repeats 
    plot.max <- set.max.repeat.size 
    
    #print(sequences$fasta.name[sequences$file.name == unique(sequences$file.name)[i]] )
    
    #print(nchar(fasta.sequence[sequences$file.name == unique(sequences$file.name)[i]]) )
    
    circos_plotter(outputs.directory, 
                   assemblyName = unique(sequences$file.name)[i], 
                   chr_names = sequences$fasta.name[sequences$file.name == unique(sequences$file.name)[i]] , 
                   chr_lengths = nchar(fasta.sequence[sequences$file.name == unique(sequences$file.name)[i]]) , 
                   plot.min, 
                   plot.max)
    
    if(simpleplot) #old plotting, can be turned on additionally
    {
      if(draw.scaffold.repeat.plots(temp.folder = execution.path, 
                                    assemblyName = strsplit(fasta.list[i], split = "/")[[1]][length(strsplit(fasta.list[i], split = "/")[[1]])], 
                                    fastaDirectory = fasta.list[i], 
                                    only.Chr1_5 = FALSE, 
                                    single.pngs = TRUE) != 0)
      {
        print("plotting likely failed")
      }
    }
    
    print(paste("finished", sep = ""))
  }
  print("Plotting finished)")
  quit()
  print("you can't see me")
}


date.now = str_replace_all(paste(strsplit(as.character(Sys.time()), split = " ")[[1]][1:2], collapse = ""), regex("[:,-]") , "")

if(!is.null(random.seed))
{
  set.seed(seed = random.seed)
} else
{
  set.seed(seed = NULL)
  random.seed = .Random.seed[3]
  set.seed(seed = random.seed)
}

write(paste("random seed is: ", random.seed,  sep = ""), 
      file = paste("TRASH_", date.now,".out", sep = ""), append = TRUE)










#identify repeats per chromosome
foreach(i = 1 : length(fasta.sequence)) %dopar% {
  #for(i in 1 : length(fasta.sequence)) {  #TODO go back to par
  print(Sys.time())
  print(paste("Working on sequence ", i , sep = ""))
  if(!skip.repetitive.regions)
  {
    repetitive.regions = NULL
    repetitive.regions = Repeat.Identifier(DNA.sequence = fasta.sequence[i], assemblyName = sequences$file.name[i], fasta.name = sequences$fasta.name[i], 
                                           kmer = set.kmer, window = window.size, threshold = set.threshold, mask.small.regions = filter.small.regions, mask.small.repeats = filter.small.repeats,
                                           max.repeat.size = set.max.repeat.size, LIMIT.REPEATS.TO.ALIGN = LIMIT.REPEATS.TO.ALIGN,
                                           tests = 4, temp.folder = execution.path, sequence.template = sequence.templates, mafft.bat.file = paste(installation.path, "/src/mafft-linux64/mafft.bat", sep = ""))
    if(typeof(repetitive.regions) == "list")
    {
      if(nrow(repetitive.regions) != 0)
      {
        write.csv(x = repetitive.regions, file = paste(execution.path, "/", sequences$file.name[i], "_out/Regions_", sequences$file.name[i], "_", sequences$fasta.name[i], ".csv", sep = ""), row.names = FALSE)
      }
    }
  }
  gc()
  
  print(paste("Finished working on sequence ", i , sep = ""))
} 









#save repeats per genome
for(i in 1 : length(fasta.list))
{
  print(Sys.time())
  print("saving repeats")
  setwd(paste(execution.path, "/", sep = ""))
  repeat.files = system(paste("find ./", unique(sequences$file.name)[i], "_out", " -name \"Regions*\"", sep = ""), intern = TRUE)
  regions = NULL
  for(j in 1 : length(repeat.files))
  {
    region.read = read.csv(file = paste(repeat.files[j], sep = ""))
    if(nrow(region.read) > 1)
    {
      regions = rbind(regions, region.read)
    }
  }
  regions = regions[,-c(1)]
  write.csv(x = regions, file = paste(execution.path, "/Summary.of.repetitive.regions.", unique(sequences$file.name)[i], ".csv", sep = ""), quote = FALSE)
  
  
  extract.all.repeats(temp.folder = execution.path, assemblyName = strsplit(fasta.list[i], split = "/")[[1]][length(strsplit(fasta.list[i], split = "/")[[1]])])
  
  export.gff(temp.folder = execution.path, 
             assemblyName = unique(sequences$file.name)[i])
  
  if(delete_temp_output)
  {
    if(file.exists(paste(execution.path, "/", sequences$file.name[i], "_out", sep = "")))
    {
      unlink(x = paste(execution.path, "/", sequences$file.name[i], "_out", sep = ""), recursive = TRUE, force = FALSE)
    }
  }
  print(paste("finished saving", sep = ""))
}


for(i in 1 : length(fasta.list))
{
  print(paste("plotting", sep = ""))
  
  outputs.directory <- paste(execution.path, sep = "")
  plot.min <- filter.small.repeats 
  plot.max <- set.max.repeat.size 
  
  #print(sequences$fasta.name[sequences$file.name == unique(sequences$file.name)[i]] )
  
  #print(nchar(fasta.sequence[sequences$file.name == unique(sequences$file.name)[i]]) )
  
  circos_plotter(outputs.directory, 
                 assemblyName = unique(sequences$file.name)[i], 
                 chr_names = sequences$fasta.name[sequences$file.name == unique(sequences$file.name)[i]] , 
                 chr_lengths = nchar(fasta.sequence[sequences$file.name == unique(sequences$file.name)[i]]) , 
                 plot.min, 
                 plot.max)
  if(simpleplot) #old plotting, can be turned on additionally
  {
    if(draw.scaffold.repeat.plots(temp.folder = execution.path, 
                                  assemblyName = strsplit(fasta.list[i], split = "/")[[1]][length(strsplit(fasta.list[i], split = "/")[[1]])], 
                                  fastaDirectory = fasta.list[i], 
                                  only.Chr1_5 = FALSE, 
                                  single.pngs = TRUE) != 0)
    {
      print("plotting likely failed")
    }
  }
  
  print(paste("finished plotting", sep = ""))
}



for(i in 1 : length(fasta.list))
{
  print(paste("edit distance", sep = ""))
  
  #do edit distance per genome
  ### this also initiates the repetitiveness column in the repeats data frame
  calc.edit.distance(temp.folder = execution.path, 
                     assemblyName = sequences$file.name[i],
                     fasta.name = sequences$fasta.name[i],
                     mafft.bat.file = paste(installation.path, "/src/mafft-linux64/mafft.bat", sep = ""))
  
  
}

for(i in 1 : length(fasta.sequence))
{
  #plot edit per chromosome
  plot.edit(temp.folder = execution.path, 
                      assemblyName = sequences$file.name[i], 
                      chr.name = sequences$fasta.name[i],
                      hor.c.script.path = hor.c.script.path,
                      class.name = class.name.for.HOR)
}


if(class.name.for.HOR != "")
{
  
  #do HORs per chromosome
  foreach(i = 1 : length(fasta.sequence)) %dopar% {
    HOR.wrapper(threshold = max.divergence.value, 
                cutoff = min.hor.value, 
                temp.folder = execution.path, 
                assemblyName = sequences$file.name[i], 
                chr.name = sequences$fasta.name[i],
                mafft.bat.file = paste(installation.path, "/src/mafft-linux64/mafft.bat", sep = ""),
                hor.c.script.path = hor.c.script.path,
                class.name = class.name.for.HOR)
    gc()
  } 
  
 

  for(i in 1 : length(fasta.sequence))
  {
    #do repetitiveness per chromosome
    calc.repetitiveness(temp.folder = execution.path, 
                        assemblyName = sequences$file.name[i], 
                        chr.name = sequences$fasta.name[i],
                        hor.c.script.path = hor.c.script.path,
                        class.name = class.name.for.HOR)
    #do plot.repetitiveness per chromosome
    plot.repetitiveness(temp.folder = execution.path, 
                        assemblyName = sequences$file.name[i], 
                        chr.name = sequences$fasta.name[i],
                        hor.c.script.path = hor.c.script.path,
                        class.name = class.name.for.HOR)
  }
}


print("TRASH finished, exiting")

gc()





