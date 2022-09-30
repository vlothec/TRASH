extract.chromosome.name = function(regname_string, assemblyName)
{
  return(strsplit(regname_string, split = paste(assemblyName, "_", sep = ""))[[1]][2])
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
  
  files = list.files(path = paste(paste(assemblyName, "_out", sep = ""), "/temp2", sep = ""), pattern = "Repeats_*")
  if(length(files) == 0)
  {
    return(1)
  }
  print("saving regions no:")
  print(length(files))
  for(i in 1 : length(files))
  {
    region = read.csv(file = paste(paste(assemblyName, "_out", sep = ""), "/temp2/", files[i], sep = ""), header = TRUE)
    print("save region:")
    print(i)
    print(nrow(region))
    if(nrow(region) > 0)
    {
      if(Sys.info()['sysname'] == "Linux")
      {
        region = region[,-c(1)]
      } 
      region.coor = strsplit(files[i], split = "_")[[1]]
      region.start = as.numeric(region.coor[length(region.coor)-1])
      region.end = as.numeric(strsplit(region.coor[length(region.coor)], split = "\\.")[[1]][1])
      region$start = region$start + region.start
      region$end = region$end + region.start
      
      names.table = names(region)
      
      write.table(x = region, sep = ",", file = paste(temp.folder, "/temp.all.repeats.from.", assemblyName, ".csv", sep = ""), row.names = FALSE, append = TRUE, col.names = FALSE)
      
      #repeats.all = rbind(repeats.all, region)
    }
  }
  repeats.all = read.csv(file = paste(temp.folder, "/temp.all.repeats.from.", assemblyName, ".csv", sep = ""), header = FALSE)
  
  names(repeats.all) = names.table
  
  repeats.all = repeats.all[order(repeats.all$start),]
  repeats.all = repeats.all[order(repeats.all$region.name),]
  
  repeats.all$seq.name = ""
  
  repeats.all$seq.name = unlist(lapply(X = repeats.all$region.name, extract.chromosome.name, assemblyName))
  
  write.csv(x = repeats.all, file = paste(temp.folder, "/all.repeats.from.", assemblyName, ".csv", sep = ""), row.names = FALSE)
  return(0)
}
