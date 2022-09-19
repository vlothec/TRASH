

plot.edit = function(temp.folder = "", #execution.path
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
  
  write(paste("plotting edit ", sep = ""), 
        file = paste(paste(assemblyName, "_out", sep = ""), "/", chr.name, ".out.txt", sep = ""), append = TRUE)
  
  repeats = read.csv(file = paste(execution.path, "/all.repeats.from.", assemblyName, ".csv", sep = ""))
  
  repeats.row.initial = nrow(repeats)
  
  classes = unique(repeats$class)
  
  print(paste("Available classes: ", classes, sep = ""))
  
  
  if(length(classes) == 0)
  {
    print(paste("no classes assigned to the repeats to plot edit", sep = ""))
    return(1)
  }
  
  repeats$class[is.na(repeats$class)] = ""
  
  repeats$seq.name[is.na(repeats$seq.name)] = ""
  
  pdf(file = paste("./plots/plot_edit_", chr.name, ".pdf", sep = ""), width = 14, height = 7)
  plot(NULL, xlim = c(min(repeats$start), max(repeats$end)), ylim = c(0,max(repeats$edit.distance[!is.na(repeats$edit.distance)])), xlab = "coordinates, bp", ylab = "edit distance")
  
  for(i in 1 : length(classes))
  {
    chromosome.repeats.IDs = which(repeats$class == classes[i] & repeats$seq.name == chr.name)
    
    
    if(length(chromosome.repeats.IDs) > 0)
    {
      
      repeats.temp = repeats[chromosome.repeats.IDs,]
      
      #plot
      lines(repeats.temp$start, movav(repeats.temp$edit.distance, nrow(repeats.temp)/50))
    }
  }
  
  
  dev.off()
 
  
  return(0)
}
