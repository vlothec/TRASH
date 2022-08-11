

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
    return(0)
  }
  for(i in 1 : length(files))
  {
    region = read.csv(file = paste(paste(assemblyName, "_out", sep = ""), "/temp2/", files[i], sep = ""), header = TRUE)
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
      repeats.all = rbind(repeats.all, region)
    }
  }
  repeats.all = repeats.all[order(repeats.all$start),]
  repeats.all = repeats.all[order(repeats.all$region.name),]
  write.csv(x = repeats.all, file = paste(temp.folder, "/all.repeats.from.", assemblyName, ".csv", sep = ""), row.names = FALSE)
  return(1)
}
