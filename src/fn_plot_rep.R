

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
  
  if(length(chromosome.repeats.IDs) > 0)
  {
    
    repeats.temp = repeats[chromosome.repeats.IDs,]
    
    #plot
    
    pdf(file = paste("./plots/plot_repetitiveness_", chr.name, ".pdf", sep = ""), width = 14, height = 7)
    plot(repeats.temp$start, movav(repeats.temp$repetitiveness, nrow(repeats.temp)/50), type = "l")
    dev.off()
    
    
  }
  
 
  
  return(0)
}
