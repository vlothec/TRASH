#!/usr/bin/env Rscript
print("start")
if(Sys.info()['sysname'] == "Linux")
{
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
} else
{
  thisFile = function() #needs to be run from the TRASH\src directory if on windows
  {
    return(strsplit(getwd(), split = "//src")[[1]][1])
  }
}



######### 
installation.path = thisFile()
installation.path = strsplit(installation.path, split = "/")[[1]]
if(Sys.info()['sysname'] == "Linux")
{
  installation.path = paste(installation.path[1:(length(installation.path) - 2)], collapse = "/")
} else
{
  installation.path = paste(installation.path[1:(length(installation.path) - 1)], collapse = "/")
}
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
  max.divergence.value = 25 # max variants allowed for HOR identification
  min.hor.value = 3 #HORs of shorter sizes than this in monomers will not be saved 
  skip.short.fasta.sequences = 0
  set.kmer = 10 # kmer size used for initial identification of repetitive regions
  set.threshold = 5 # window repetitiveness score (0-100) threshold
  set.max.repeat.size = 850 # max size of repeats to be identified
  filter.small.regions = 3000 # repetitive windows smaller than this size will be removed (helps getting rid of regions with short duplications)
  filter.small.repeats = 5 # repetitive windows where dominant kmer distance is lower than this value will be removed (for example AT dinucleotide repeats)
  window.size = 1000 # how far apart kmers can be in the initial search for exact matches. No repeats larger than this will be identified
  MAX.CHROMOSOMES.TODO = 48 
  hor.only = FALSE
  class.name.for.HOR = ""
  delete_temp_output = FALSE
  LIMIT.REPEATS.TO.ALIGN = 35600 #in base pairs
  simpleplot = TRUE
  circosplot = FALSE
  random.seed = NULL
  N.max.div = 100 #fn_new_distance threshold score above which will look for divisions, the lower, the more loose
  try.until = 12 #fn_new_distance max number of N divisions, the higher the longer repeats can be split
  smooth.percent = 2 #fn_new_distance smoothing factor for finding the histogram peaks, the higher the wider
  plot.N = F #whether N dists should be plotted, TODO add flag
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
      else if(arguments[i] == "--circosplot")
      {
        circosplot = TRUE
      }
      else if(arguments[i] == "--par")
      {
        cores.no = as.numeric(arguments[i + 1])
      }
      else if(arguments[i] == "--N.max.div")
      {
        N.max.div = as.numeric(arguments[i + 1])
      }
      else if(arguments[i] == "--max.N.split")
      {
        try.until = as.numeric(arguments[i + 1])
      }
      else if(arguments[i] == "--smooth.percent")
      {
        smooth.percent = as.numeric(arguments[i + 1])
      }
      if(grepl("\\.(fa|fna|fasta)$", arguments[i]))
      {
        if(file.exists(paste(execution.path, arguments[i], sep = "/")))
        {
          if(!dir.exists(paste(execution.path, arguments[i], sep = "/")))
          {
            if(startsWith(readLines(paste(execution.path, arguments[i], sep = "/"), n = 1), ">"))
            {
              fasta.list = c(fasta.list, paste(execution.path, arguments[i], sep = "/"))
            }
          }
        } else if(file.exists(arguments[i]))
        {
          if(!dir.exists(arguments[i]))
          { 
            if(startsWith(readLines(arguments[i], n = 1), ">"))
            {
              fasta.list = c(fasta.list, arguments[i])
            }
          }
        } else
        {
          warning(paste("Cannot find the", arguments[i], "fasta file, proceeding regardless", sep = " "), call. = FALSE)
        }
      }
    }
  }
  if(Sys.info()['sysname'] != "Linux")
  {
    print(Sys.info()['sysname'])
    print("Only Linux OS tested")
  } 
}

if(length(fasta.list) == 0) stop("No fasta files found, exiting", call. = FALSE)

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
suppressPackageStartupMessages(library(circlize, quietly = TRUE))
suppressPackageStartupMessages(library(pwalign, quietly = TRUE))

if(Sys.info()['sysname'] == "Linux")
{
  mafft.local = "/src/mafft-linux64/mafft.bat"
} else
{
  mafft.local = "/src/mafft-win/mafft.bat"
}



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
        if(Sys.info()['sysname'] == "Windows")
        {
          sequences = rbind(sequences, data.frame(file.name = strsplit(fasta.list[i], split = "\\\\")[[1]][length(strsplit(fasta.list[i], split = "\\\\")[[1]])], fasta.name = names(fasta[j])))
          
        }else
        {
          sequences = rbind(sequences, data.frame(file.name = strsplit(fasta.list[i], split = "/")[[1]][length(strsplit(fasta.list[i], split = "/")[[1]])], fasta.name = names(fasta[j])))
          
        }
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
  
  # for(i in 1 : length(fasta.list))
  # {
  #   print(paste("edit distance", sep = ""))
  #   
  #   #do edit distance per genome
  #   ### this also initiates the repetitiveness column in the repeats data frame
  #   calc.edit.distance(temp.folder = execution.path, 
  #                      assemblyName = sequences$file.name[i],
  #                      fasta.name = sequences$fasta.name[i],
  #                      LIMIT.REPEATS.TO.ALIGN = LIMIT.REPEATS.TO.ALIGN,
  #                      mafft.bat.file = paste(installation.path, mafft.local, sep = ""))
  #   
  #   
  # }
  # gc()
  # for(i in 1 : length(fasta.sequence))
  # {
  #   #plot edit per chromosome
  #   plot.edit(temp.folder = execution.path, 
  #             assemblyName = sequences$file.name[i], 
  #             chr.name = sequences$fasta.name[i],
  #             hor.c.script.path = hor.c.script.path,
  #             class.name = class.name.for.HOR)
  # }
  # gc()
  
  
  print("calculating HORs")
  
  #do HORs per chromosome
  if(Sys.info()['sysname'] == "Windows")
  {
    
    for(i in 1 : length(fasta.sequence)) {
      HOR.wrapper(threshold = max.divergence.value, 
                  cutoff = min.hor.value,
                  temp.folder = execution.path, 
                  assemblyName = sequences$file.name[i], 
                  chr.name = sequences$fasta.name[i],
                  mafft.bat.file = paste(installation.path, mafft.local, sep = ""),
                  hor.c.script.path = hor.c.script.path,
                  class.name = class.name.for.HOR)
      gc()
    } 
  } else
  {
    foreach(i = 1 : length(fasta.sequence)) %dopar% {  
      HOR.wrapper(threshold = max.divergence.value, 
                  cutoff = min.hor.value,
                  temp.folder = execution.path, 
                  assemblyName = sequences$file.name[i], 
                  chr.name = sequences$fasta.name[i],
                  mafft.bat.file = paste(installation.path, mafft.local, sep = ""),
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
    plot.min <- filter.small.repeats #TODO add a flag to change
    plot.max <- set.max.repeat.size  #TODO add a flag to change
    #TODO add a flag to change threshold
    
    #print(sequences$fasta.name[sequences$file.name == unique(sequences$file.name)[i]] )
    
    #print(nchar(fasta.sequence[sequences$file.name == unique(sequences$file.name)[i]]) )
    if(circosplot) {
      circos_plotter(outputs.directory, 
                     assemblyName = unique(sequences$file.name)[i], 
                     chr_names = sequences$fasta.name[sequences$file.name == unique(sequences$file.name)[i]] , 
                     chr_lengths = nchar(fasta.sequence[sequences$file.name == unique(sequences$file.name)[i]]) , 
                     plot.min, 
                     plot.max)
    }
    
    
    if(simpleplot) #old plotting, can be turned on additionally
    {
      if(draw.scaffold.repeat.plots(temp.folder = execution.path, 
                                    assemblyName = strsplit(fasta.list[i], split = "/")[[1]][length(strsplit(fasta.list[i], split = "/")[[1]])], 
                                    fastaDirectory = fasta.list[i], 
                                    single.pngs = TRUE, y.max = set.max.repeat.size) != 0)
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



print("Sequences to analyse:")

print(sequences)

#identify repeats per chromosome
if(Sys.info()['sysname'] == "Windows")
{
  for(i in 1 : length(fasta.sequence)) {  
    print(Sys.time())
    print(paste("Working on sequence ", i , sep = ""))
    if(!skip.repetitive.regions)
    {
      repetitive.regions = NULL
      repetitive.regions = Repeat.Identifier(DNA.sequence = fasta.sequence[i], assemblyName = sequences$file.name[i], fasta.name = sequences$fasta.name[i], 
                                             kmer = set.kmer, window = window.size, threshold = set.threshold, mask.small.regions = filter.small.regions, mask.small.repeats = filter.small.repeats,
                                             max.repeat.size = set.max.repeat.size, LIMIT.REPEATS.TO.ALIGN = LIMIT.REPEATS.TO.ALIGN,
                                             tests = 6, temp.folder = execution.path, sequence.template = sequence.templates, mafft.bat.file = paste(installation.path, mafft.local, sep = ""),
                                             N.max.div , try.until, smooth.percent, plot.N = plot.N)
      if(typeof(repetitive.regions) == "list")
      {
        if(length(repetitive.regions) != 0)
        {
          write.csv(x = repetitive.regions, file = paste(execution.path, "/", sequences$file.name[i], "_out/Regions_", sequences$file.name[i], "_", sequences$fasta.name[i], ".csv", sep = ""), row.names = FALSE)
        }
      }
    }
  }
} else
{
  foreach(i = 1 : length(fasta.sequence)) %dopar% {
    print(Sys.time())
    print(paste("Working on sequence ", i , sep = ""))
    if(!skip.repetitive.regions)
    {
      repetitive.regions = NULL
      repetitive.regions = Repeat.Identifier(DNA.sequence = fasta.sequence[i], assemblyName = sequences$file.name[i], fasta.name = sequences$fasta.name[i], 
                                             kmer = set.kmer, window = window.size, threshold = set.threshold, mask.small.regions = filter.small.regions, mask.small.repeats = filter.small.repeats,
                                             max.repeat.size = set.max.repeat.size, LIMIT.REPEATS.TO.ALIGN = LIMIT.REPEATS.TO.ALIGN,
                                             tests = 6, temp.folder = execution.path, sequence.template = sequence.templates, mafft.bat.file = paste(installation.path, mafft.local, sep = ""),
                                             N.max.div , try.until, smooth.percent, plot.N = plot.N)
      if(typeof(repetitive.regions) == "list")
      {
        if(length(repetitive.regions) != 0)
        {
          write.csv(x = repetitive.regions, file = paste(execution.path, "/", sequences$file.name[i], "_out/Regions_", sequences$file.name[i], "_", sequences$fasta.name[i], ".csv", sep = ""), row.names = FALSE)
        }
      }
    }
  }
} 
gc()

print(paste("Finished working on sequence ", i , sep = ""))









#save repeats per genome
for(i in 1 : length(fasta.list))
{
  print(Sys.time())
  print("saving repeats")
  setwd(paste(execution.path, "/", sep = ""))
  print(execution.path)
  repeat.files = list.files(path = paste("./", unique(sequences$file.name)[i], "_out", sep = ""), pattern = "Regions*", full.names = T)
  if(length(repeat.files) > 0)
  {
    print(repeat.files)
    #repeat.files = system(paste("find ./", unique(sequences$file.name)[i], "_out", " -name \"Regions*\"", sep = ""), intern = TRUE)
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
    
    
    extract.all.repeats(temp.folder = execution.path, assemblyName = unique(sequences$file.name)[i])
    
    gc()
    
    export.gff(temp.folder = execution.path, 
               assemblyName = unique(sequences$file.name)[i])
    
    gc()
    
    if(delete_temp_output)
    {
      if(file.exists(paste(execution.path, "/", sequences$file.name[i], "_out", sep = "")))
      {
        unlink(x = paste(execution.path, "/", sequences$file.name[i], "_out", sep = ""), recursive = TRUE, force = FALSE)
      }
    }
    print(paste("finished saving", sep = ""))
  }
}



for(i in 1 : length(fasta.list))
{
  print(paste("edit distance", sep = ""))
  
  #do edit distance per genome
  ### this also initiates the repetitiveness column in the repeats data frame
  calc.edit.distance(temp.folder = execution.path, 
                     assemblyName = sequences$file.name[i],
                     fasta.name = sequences$fasta.name[i],
                     LIMIT.REPEATS.TO.ALIGN = LIMIT.REPEATS.TO.ALIGN,
                     mafft.bat.file = paste(installation.path, mafft.local, sep = ""))
  
  
}
gc()
for(i in 1 : length(fasta.sequence))
{
  #plot edit per chromosome
  plot.edit(temp.folder = execution.path, 
            assemblyName = sequences$file.name[i], 
            chr.name = sequences$fasta.name[i],
            hor.c.script.path = hor.c.script.path,
            class.name = class.name.for.HOR)
}
gc()

if(class.name.for.HOR != "")
{
  
  #do HORs per chromosome
  if(Sys.info()['sysname'] == "Windows")
  {
    for(i in 1 : length(fasta.sequence)) {
      HOR.wrapper(threshold = max.divergence.value, 
                  cutoff = min.hor.value, 
                  temp.folder = execution.path, 
                  assemblyName = sequences$file.name[i], 
                  chr.name = sequences$fasta.name[i],
                  mafft.bat.file = paste(installation.path, mafft.local, sep = ""),
                  hor.c.script.path = hor.c.script.path,
                  class.name = class.name.for.HOR)
    }
  } else
  {
    foreach(i = 1 : length(fasta.sequence)) %dopar% {
      HOR.wrapper(threshold = max.divergence.value, 
                  cutoff = min.hor.value, 
                  temp.folder = execution.path, 
                  assemblyName = sequences$file.name[i], 
                  chr.name = sequences$fasta.name[i],
                  mafft.bat.file = paste(installation.path, mafft.local, sep = ""),
                  hor.c.script.path = hor.c.script.path,
                  class.name = class.name.for.HOR)
    } 
  }
  gc()
  
  
  
  
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
                                  single.pngs = TRUE, y.max = set.max.repeat.size) != 0)
    {
      print("plotting likely failed")
    }
  }
  
  print(paste("finished plotting", sep = ""))
}


print("TRASH finished, exiting")

gc()





