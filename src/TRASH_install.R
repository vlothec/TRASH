#!/usr/bin/env Rscript
cmd.Args = commandArgs(trailingOnly = TRUE)
change.lib.paths = TRUE
if(length(cmd.Args) > 0)
{
  if(cmd.Args[1] == "--def")
  {
    change.lib.paths = FALSE
  } 
}

if(as.numeric(strsplit(strsplit(R.version.string, ' ')[[1]][3], "[.]")[[1]][1]) < 4)
{
  stop("R version found is lower than 4")
}

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

print("inst path:")
print(installation.path)
lib.path = paste(installation.path, "/libs", sep = "")

print("lib.path:")
print(lib.path)

if(!dir.exists(lib.path))
{
  dir.create(lib.path)
}

if(change.lib.paths)
{
  .libPaths(lib.path)
}

print("loading remotes")
if(!require("remotes", quietly = TRUE, warn.conflicts = FALSE, lib.loc = .libPaths()[1]))
{
  install.packages("remotes", repos = "https://cloud.r-project.org", lib = .libPaths()[1])
  library("remotes", quietly = TRUE, warn.conflicts = FALSE, lib.loc = .libPaths()[1])
}

print("loading stringr")
if(!require("stringr", quietly = TRUE, warn.conflicts = FALSE,  lib.loc = .libPaths()[1]))
{
  if(Sys.info()['sysname'] == "Linux")
  {
    install_version("stringr", version = "1.4.0", repos = "http://cran.us.r-project.org")
  } else
  {
    install.packages("stringr", repos = "http://cran.us.r-project.org", lib = .libPaths()[1])
  }
  library("stringr", quietly = TRUE, warn.conflicts = FALSE, lib.loc = .libPaths()[1])
  if(packageVersion("stringr") != "1.4.0")
  {
    print(paste("\"stringr\" library version is different than recommended (1.4.0). Consider installing TRASH in a new folder (see manual)"))
  }
}

print("loading stringdist")
if(!require("stringdist", quietly = TRUE, warn.conflicts = FALSE, lib.loc = .libPaths()[1]))
{
  if(Sys.info()['sysname'] == "Linux")
  {
    install_version("stringdist", version = "0.9.8", repos = "http://cran.us.r-project.org")
  } else
  {
    install.packages("stringdist", repos = "http://cran.us.r-project.org", lib = .libPaths()[1])
  }
  
  library("stringdist", quietly = TRUE, warn.conflicts = FALSE, lib.loc = .libPaths()[1])
  if(packageVersion("stringdist") != "0.9.8")
  {
    print(paste("\"stringdist\" library version is different than recommended (0.9.8). Consider installing TRASH in a new folder (see manual)"))
  }
}

print("loading base")
if(!require("base"))
{
  if(Sys.info()['sysname'] == "Linux")
  {
    install_version("base", version = "4.0.3", repos = "http://cran.us.r-project.org")
  } else
  {
    install.packages("base", repos = "http://cran.us.r-project.org", lib = .libPaths()[1])
  }
  
  library("base", quietly = TRUE, warn.conflicts = FALSE, lib.loc = .libPaths()[1])
  if(packageVersion("base") != "4.0.3")
  {
    print(paste("\"base\" library version is different than recommended (4.0.3). Consider installing TRASH in a new folder (see manual)"))
  }
}

print("loading BiocManager")
if(!require("BiocManager", quietly = TRUE, warn.conflicts = FALSE, lib.loc = .libPaths()[1]))
{
  if(Sys.info()['sysname'] == "Linux")
  {
    install.packages("BiocManager", version = "1.30.16", repos = "http://cran.us.r-project.org")
  } else
  {
    install.packages("BiocManager", repos = "http://cran.us.r-project.org", lib = .libPaths()[1])
  }
  
  if(packageVersion("BiocManager") != "1.30.16")
  {
    print(paste("\"BiocManager\" library version is different than recommended (1.30.16). Consider installing TRASH in a new folder (see manual)"))
  }
  BiocManager::install()
}

print("loading Biostrings")
if(!require("Biostrings", quietly = TRUE, warn.conflicts = FALSE, lib.loc = .libPaths()[1]))
{
  BiocManager::install("Biostrings")
  
  library("Biostrings", quietly = TRUE, warn.conflicts = FALSE, lib.loc = .libPaths()[1])
}

print("loading circlize")
if(!require("circlize", quietly = TRUE, warn.conflicts = FALSE, lib.loc = .libPaths()[1]))
{
  if(Sys.info()['sysname'] == "Linux")
  {
    install_version("circlize", version = "0.4.15", repos = "http://cran.us.r-project.org")
  } else
  {
    install.packages("circlize", repos = "http://cran.us.r-project.org", lib = .libPaths()[1])
  }
  
  library("circlize", quietly = TRUE, warn.conflicts = FALSE, lib.loc = .libPaths()[1])
  if(packageVersion("circlize") != "0.4.15")
  {
    print(paste("\"circlize\" library version is different than recommended (0.4.15). Consider installing TRASH in a new folder (see manual)"))
  }
}

print("loading seqinr")
if(!require("seqinr", quietly = TRUE, warn.conflicts = FALSE, lib.loc = .libPaths()[1]))
{
  if(Sys.info()['sysname'] == "Linux")
  {
    install_version("seqinr", version = "4.2.8", repos = "http://cran.us.r-project.org")
  } else
  {
    install.packages("seqinr", repos = "http://cran.us.r-project.org", lib = .libPaths()[1])
  }
  
  library("seqinr", quietly = TRUE, warn.conflicts = FALSE, lib.loc = .libPaths()[1])
  if(packageVersion("seqinr") != "4.2.8")
  {
    print(paste("\"seqinr\" library version is different than recommended (4.2.8). Consider installing TRASH in a new folder (see manual)"))
  }
}

print("loading doParallel")
if(!require("doParallel", quietly = TRUE, warn.conflicts = FALSE, lib.loc = .libPaths()[1]))
{
  if(Sys.info()['sysname'] == "Linux")
  {
    install_version("doParallel", version = "1.0.17", repos = "http://cran.us.r-project.org")
  } else
  {
    install.packages("doParallel", repos = "http://cran.us.r-project.org", lib = .libPaths()[1])
  }
  
  library("doParallel", quietly = TRUE, warn.conflicts = FALSE, lib.loc = .libPaths()[1])
  if(packageVersion("doParallel") != "1.0.17")
  {
    print(paste("\"doParallel\" library version is different than recommended (1.0.17). Consider installing TRASH in a new folder (see manual)"))
  }
}

src.files = list.files(path = paste(installation.path, "/src", sep = ""), pattern = "fn_", full.names = TRUE)

for(i in 1 : length(src.files))
{
  print(paste("attaching function from: ", src.files[i], sep = ""))
  source(src.files[i])
}

print("Checking mafft installation")
if(Sys.info()['sysname'] == "Linux")
{
  if(file.exists(paste(installation.path, "/src/mafft-linux64/mafft.bat", sep = "")))
  {
    print("mafft executable found")
  } else if(file.exists(paste(installation.path, "/src/mafft-7.490-linux.tgz", sep = "")))
  {
    print("Unpacking mafft")
    system(paste("tar -xzvf ", installation.path, "/src/mafft-7.490-linux.tgz", sep = ""), intern = TRUE)
  } else 
  {
    system(paste("wget https://mafft.cbrc.jp/alignment/software/mafft-7.490-linux.tgz", sep = ""))
    print("Unpacking mafft")
    system(paste("tar -xzvf ", installation.path, "/mafft-7.490-linux.tgz -C ", installation.path, "/src", sep = ""), intern = TRUE)
  }
  
  if(!file.exists(paste(installation.path, "/src/mafft-linux64/mafft.bat", sep = "")))
  {
    warning("neither mafft executable nor package found \nmake sure mafft is downloaded in the /src directory \nhttps://mafft.cbrc.jp/alignment/software/mafft-7.490-linux.tgz\nor unable to download")
  }
} else
{
  if(file.exists(paste(installation.path, "/src/mafft-win/mafft.bat", sep = "")))
  {
    print("mafft executable found")
  }  else 
  {
    download.file("https://mafft.cbrc.jp/alignment/software/mafft-7.511-win64-signed.zip", "mafft.zip")
    print("Unpacking mafft")
    unzip(zipfile = "mafft.zip")
  }
  
  if(!file.exists(paste(installation.path, "/src/mafft-win/mafft.bat", sep = "")))
  {
    warning("neither mafft executable nor package found \nmake sure mafft is downloaded in the /src directory \nhttps://mafft.cbrc.jp/alignment/software/mafft-7.490-linux.tgz\nor unable to download")
  }
}



print("Checking R library installation")

if(exists("draw.scaffold.repeat.plots") & 
   exists("extract.all.repeats") & 
   exists("extract_kmers") & 
   exists("Hash_And_Reverse") & 
   exists("kmer.compare") & 
   exists("revCompString")){
  print("All functions available and libraries installed, exiting")
}











