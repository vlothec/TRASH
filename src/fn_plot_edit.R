

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
  
  if(!file.exists(paste(execution.path, "/all.repeats.from.", assemblyName, ".csv", sep = "")))
  {
    write(paste("No repeats in the file, edits not plotted", sep = ""), 
          file = paste(paste(assemblyName, "_out", sep = ""), "/", chr.name, ".out.txt", sep = ""), append = TRUE)
    return(1)
  }
  
  
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
  
  colour.roll = rep(pal <- c("#000000","#004949","#009292","#ff6db6","#ffb6db",
                             "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
                             "#920000","#924900","#db6d00","#24ff24","#ffff6d"), 
                    ceiling(length(classes)/15))
  
  pdf(file = paste("./plots/plot_edit_", chr.name, ".pdf", sep = ""), width = 14, height = 7)
  plot(NULL, 
       xlim = c(min(repeats$start), max(repeats$end)), 
       #ylim = c(0,max(repeats$edit.distance[!is.na(repeats$edit.distance)])), 
       ylim = c(0,100), #0 to 100 because edit distance in normalized
       xlab = "coordinates, bp", 
       ylab = "edit distance normalised to max possible",
       main = paste("edit distances ", chr.name, sep = ""))
  for(i in 1 : length(classes))
  {
    chromosome.repeats.IDs = which(repeats$class == classes[i] & repeats$seq.name == chr.name)
    
    
    if(length(chromosome.repeats.IDs) > 0)
    {
      
      repeats.temp = repeats[chromosome.repeats.IDs,]
      
      windows = seq(min(repeats.temp$start), max(repeats.temp$start), by = 10000)
      
      windows.edit = NULL
      
      for(k in 1:length(windows))
      {
        windows.edit = c(windows.edit, mean(repeats.temp$edit.distance[which(repeats.temp$start > windows[k] & repeats.temp$start <= windows[k + 1])]))
      }
      
      #windows.edit.norm = 100 *  windows.edit / quantile(windows.edit, probs = 0.999, na.rm = T)
      windows.edit.norm = 100 *  windows.edit / ceiling(mean(repeats.temp$width[!is.na(repeats.temp$width)]))
      
      #plot
      lines(windows, windows.edit.norm, type = "o", col = colour.roll[i])
    }
  }
  
  
  dev.off()
 
  
  return(0)
}
