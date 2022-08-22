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

installation.path = thisFile()
installation.path = strsplit(installation.path, split = "/")[[1]]
installation.path = paste(installation.path[1:(length(installation.path) - 2)], collapse = "/")
execution.path = getwd()

print("inst path:")
print(installation.path)
lib.path = paste(installation.path, "/libs", sep = "")

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
if(!require("stringr", quietly = TRUE, warn.conflicts = FALSE))
{
  install_version("stringr", version = "1.4.0", repos = "http://cran.us.r-project.org")
  library("stringr", quietly = TRUE, warn.conflicts = FALSE)
  if(packageVersion("stringr") != "1.4.0")
  {
    print(paste("\"stringr\" library version is different than recommended (1.4.0). Consider installing TRASH in a new folder (see manual)"))
  }
}

print("loading base")
if(!require("base"))
{
  install_version("base", version = "4.0.3", repos = "http://cran.us.r-project.org")
  library("base", quietly = TRUE, warn.conflicts = FALSE)
  if(packageVersion("base") != "4.0.3")
  {
    print(paste("\"base\" library version is different than recommended (4.0.3). Consider installing TRASH in a new folder (see manual)"))
  }
}

print("loading BiocManager")
if(!require("BiocManager", quietly = TRUE, warn.conflicts = FALSE))
{
  install.packages("BiocManager")
  #install_version("BiocManager", version = "1.30.16", repos = "http://cran.us.r-project.org")
  if(packageVersion("BiocManager") != "1.30.16")
  {
    print(paste("\"BiocManager\" library version is different than recommended (1.30.16). Consider installing TRASH in a new folder (see manual)"))
  }
  BiocManager::install()
}

print("loading Biostrings")
if(!require("Biostrings", quietly = TRUE, warn.conflicts = FALSE))
{
  BiocManager::install("Biostrings")
  library("Biostrings", quietly = TRUE, warn.conflicts = FALSE)
}

print("loading seqinr")
if(!require("seqinr", quietly = TRUE, warn.conflicts = FALSE))
{
  install_version("seqinr", version = "4.2.8", repos = "http://cran.us.r-project.org")
  library("seqinr", quietly = TRUE, warn.conflicts = FALSE)
  if(packageVersion("seqinr") != "4.2.8")
  {
    print(paste("\"seqinr\" library version is different than recommended (4.2.8). Consider installing TRASH in a new folder (see manual)"))
  }
}

print("loading doParallel")
if(!require("doParallel", quietly = TRUE, warn.conflicts = FALSE))
{
  install_version("doParallel", version = "1.0.17", repos = "http://cran.us.r-project.org")
  library("doParallel", quietly = TRUE, warn.conflicts = FALSE)
  if(packageVersion("doParallel") != "1.0.17")
  {
    print(paste("\"doParallel\" library version is different than recommended (1.0.17). Consider installing TRASH in a new folder (see manual)"))
  }
}

print("loading circlize")
if(!require("circlize", quietly = TRUE, warn.conflicts = FALSE))
{
  install_version("circlize", version = "0.4.15", repos = "http://cran.us.r-project.org")
  library("circlize", quietly = TRUE, warn.conflicts = FALSE)
  if(packageVersion("circlize") != "0.4.15")
  {
    print(paste("\"circlize\" library version is different than recommended (0.4.15). Consider installing TRASH in a new folder (see manual)"))
  }
}

src.files = list.files(path = paste(installation.path, "/src", sep = ""), pattern = "fn_", full.names = TRUE)

for(i in 1 : length(src.files))
{
  print(paste("attaching function from: ", src.files[i], sep = ""))
  source(src.files[i])
}

print("Checking installation")

if(exists("draw.scaffold.repeat.plots") & 
   exists("extract.all.repeats") & 
   exists("extract_kmers") & 
   exists("Hash_And_Reverse") & 
   exists("kmer.compare") & 
   exists("revCompString")){
  print("All functions available and libraries installed, exiting")
}











