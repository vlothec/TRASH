

calc.repetitiveness = function(temp.folder = "", #execution.path
                               assemblyName = "",#sequences$file.name[i], 
                               chr.name = "",#sequences$fasta.name[i],
                               hor.c.script.path = "",#hor.c.script.path,
                               class.name = "")
{
  write(paste("calculating repetitiveness", sep = ""), 
        file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
  
  repeats = read.csv(file = paste(execution.path, "/all.repeats.from.", assemblyName, ".csv", sep = ""))
  
  repeats.row.initial = nrow(repeats)
  
  HORs = read.csv(file = paste(temp.folder, "/HORs_", assemblyName, "_", chr.name, "_", class.name, ".csv", sep = ""))
  
  classes = unique(repeats$class)
  
  print(paste("Available classes: ", classes, sep = ""))
  
  
  if(length(classes) == 0)
  {
    print(paste("no classes assigned to the repeats to calculate repetitiveness", sep = ""))
    return(1)
  }
  
  repeats$class[is.na(repeats$class)] = ""
  
  repeats$seq.name[is.na(repeats$seq.name)] = ""
  
  chromosome.repeats.IDs = which(repeats$seq.name == chr.name)
  
  if(length(chromosome.repeats.IDs) > 0)
  {
    
    repeats.temp = repeats[chromosome.repeats.IDs,]
    
    repeats.temp$repetitiveness = 0
    
    for(i in 1 : nrow(HORs))
    {
      for(j in 0 : (HORs$block.size.in.units[i] - 1))
      {
        
        repeats.temp$repetitiveness[HORs$start_A[i] + j] = repeats.temp$repetitiveness[HORs$start_A[i] + j] + HORs$block.size.in.units[i]
        repeats.temp$repetitiveness[HORs$start_B[i] + j] = repeats.temp$repetitiveness[HORs$start_B[i] + j] + HORs$block.size.in.units[i]
        
      }
    }
    
  }
  
  repeats$repetitiveness[chromosome.repeats.IDs] = repeats.temp$repetitiveness
  
  #### check if the data frame is the same length as at the beginning and remove first column if necessary
  
  if(repeats.row.initial != nrow(repeats))
  {
    write(paste("diffrent number of rows in the edited repeats file, making a new file to compare", sep = ""), 
          file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
    
    write.csv(x = repeats, file = paste("temp_", repeats.file, sep = ""), row.names = FALSE)
    
  } else #### export the repeats data frame with edit distance 
  {
    write(paste("extracting repeats data frame with edit distance", sep = ""), 
          file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
    
    write.csv(x = repeats, file = paste(repeats.file, sep = ""), row.names = FALSE)
  }
  
  return(0)
}
