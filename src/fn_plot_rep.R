

plot.repetitiveness = function(temp.folder = "", #execution.path
                                        assemblyName = "",#sequences$file.name[i], 
                                        chr.name = "",#sequences$fasta.name[i],
                                        hor.c.script.path = "",#hor.c.script.path,
                                        class.name = "")
{
  if(temp.folder == "")
  {
    stop("no temp folder")
  }
  setwd(temp.folder)
  
  if(!dir.exists("./plots"))
  {
    dir.create("./plots")
  }
  
  write(paste("plotting repetitiveness", sep = ""), 
        file = paste(paste(assemblyName, "_out", sep = ""), "/", chr.name, ".out.txt", sep = ""), append = TRUE)
  
  if(!file.exists(paste(execution.path, "/all.repeats.from.", assemblyName, ".csv", sep = "")))
  {
    write(paste("No repeats in the file, repetitiveness not plotted", sep = ""), 
          file = paste(paste(assemblyName, "_out", sep = ""), "/", chr.name, ".out.txt", sep = ""), append = TRUE)
    return(1)
  }
  
  repeats = read.csv(file = paste(execution.path, "/all.repeats.from.", assemblyName, ".csv", sep = ""))
  
  repeats.row.initial = nrow(repeats)
  
  classes = unique(repeats$class)
  
  print(paste("Available classes: ", classes, sep = ""))
  
  
  if(length(classes) == 0)
  {
    print(paste("no classes assigned to the repeats to plot repetitiveness", sep = ""))
    return(1)
  }
  
  repeats$class[is.na(repeats$class)] = ""
  
  repeats$seq.name[is.na(repeats$seq.name)] = ""
  
  chromosome.repeats.IDs = which(repeats$class == class.name & repeats$seq.name == chr.name)
  
  if(length(chromosome.repeats.IDs) > 1)
  {
    
    repeats.temp = repeats[chromosome.repeats.IDs,]
    
    if(length(which(repeats.temp$repetitiveness > 0)) > 0)
    {
      print("plot repeat")
      print(chromosome.repeats.IDs)
      print(repeats.temp$repetitiveness)
      
      windows = seq(min(repeats.temp$start), max(repeats.temp$start), by = 10000)
      windows.repetitiveness = NULL
      for(k in 1:length(windows))
      {
        windows.repetitiveness = c(windows.repetitiveness, mean(repeats.temp$repetitiveness[which(repeats.temp$start > windows[k] & repeats.temp$start <= windows[k + 1])]))
      }
      print(windows.repetitiveness)
      #plot
      
      pdf(file = paste("./plots/plot_repetitiveness_", chr.name, "_", class.name, ".pdf", sep = ""), width = 14, height = 7)
      plot(windows, windows.repetitiveness, 
           type = "o", 
           main = paste("Repetitiveness ", chr.name, sep = ""),
           xlab = "coordinates, bp",
           ylab = "repetitiveness")
      dev.off()
    }
    
    
    
    
  }
  
 
  
  return(0)
}
