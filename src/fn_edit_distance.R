

calc.edit.distance = function(temp.folder = "",#execution.path
                              assemblyName = "",#sequences$file.name[i]
                              fasta.name = "",
                              mafft.bat.file = "")#paste(installation.path, "/src/mafft-linux64/mafft.bat", sep = "")
{
  write(paste("calculating edit distance", sep = ""), 
        file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
  
  repeats.file = paste(execution.path, "/all.repeats.from.", assemblyName, ".csv", sep = "")
  
  
  if(!file.exists(repeats.file))
  {
    print("cannot find repeat file, stopping edit distance calculation")
    return(1)
  }
  
  if(file.size(repeats.file) == 0)
  {
    print("input repeat file file for edit distance empty, stopping edit distance calculation")
    return(1)
  }
  
  repeats = read.csv(file = repeats.file)
  
  repeats.row.initial = nrow(repeats)
  
  if(repeats.row.initial == 0)
  {
    print("input repeat file file for edit distance has 0 rows, stopping edit distance calculation")
    return(1)
  }
  
  classes = unique(repeats$class)
  
  classes = classes[!is.na(classes)]
  
  print(paste("Available classes: ", paste(classes, collapse = ", "), sep = ""))
  
  
  if(length(classes) == 0)
  {
    print(paste("no classes assigned to the repeats to calculate edit distance, stopping edit distance calculation", sep = ""))
    return(1)
  }
  
  repeats$class[is.na(repeats$class)] = ""
  
  repeats$seq.name[is.na(repeats$seq.name)] = ""
  
  classes = unique(repeats$class)
  
  
  rm_class = NULL
  
  for(i in 1 : length(classes))
  {
    if(classes[i] == "")
    {
      rm_class = c(rm_class, i)
    }
  }
  classes = classes[-rm_class]
  
  repeats$edit.distance = NA
  repeats$repetitiveness = NA
  
  for(i in 1 : length(classes))
  {
    repeats.class.IDs = which(repeats$class == classes[i])
    
    #### log the class and how many repeats will be aligned
    
    write(paste("class: ", classes[i], " count: ", length(repeats.class.IDs), sep = ""), 
          file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
    
    if(length(repeats.class.IDs) > 1)
    {
      #### make a temp data frame with the repeats of the single class
      
      repeats.temp = repeats[repeats$class == classes[i],]
      
      #### extract and align repeats of the class
      
      file.to.align = paste(paste(assemblyName, "_out", sep = ""), "/temp1/edit.class.", i, ".fasta", sep = "")
      
      write.fasta(sequences = as.list(repeats.temp$seq), names = as.list(paste(1:length(repeats.temp$strand), "_D", repeats.temp$strand,  sep = "")), file.out = file.to.align, open = "w")
      
      if(!file.exists(file.to.align))
      {
        print("cannot find file to align for HORs, stopping HOR calculation")
        return(1)
      }
      
      if(file.size(file.to.align) == 0)
      {
        print("input fasta file for HORs empty, stopping HOR calculation")
        return(1)
      }
      
      output.alignment = paste(paste(assemblyName, "_out", sep = ""), "/temp1/aligned.edit.class.", i, ".fasta", sep = "")
      
      system(paste(mafft.bat.file, " --quiet --inputorder --kimura 1 --retree 1 ", file.to.align, " > ", output.alignment, sep = ""), intern = TRUE)
      
      if(!file.exists(output.alignment))
      {
        print("cannot find alignment file for edit distance, stopping edit distance calculation")
        return(1)
      }
      
      if(file.size(output.alignment) == 0)
      {
        print("output alignment file for edit distance empty, stopping edit distance calculation")
        return(1)
      }
      
      
      #### import the alignment
      
      alignment = read.alignment(output.alignment, format = "fasta", forceToLower = FALSE)
      
      #### generate consensus
      
      consensus.sequence = consensus(alignment, method = "majority")
      consensus.sequence = consensus.sequence[consensus.sequence != "-"]
      consensus.sequence = toupper(paste(consensus.sequence, collapse = ""))
      
      #### calculate edit distance
      
      for(j in 1 : nrow(repeats.temp))
      {
        repeats.temp$edit.distance[j] = stringdist(consensus.sequence, repeats.temp$seq[j], method = "lv")
      }
      
      #### update the repeats data frame
      
      repeats$edit.distance[repeats.class.IDs] = repeats.temp$edit.distance
      
      remove(repeats.temp)
      
    }
    #### log the end of the class edit distance calculation
    
    write(paste("class: ", classes[i], ", end of class edit distance calculation", sep = ""), 
          file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
    
  }
  
  #### log the end of all calculations
  
  write(paste("end of all classes edit distance calculation", sep = ""), 
        file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
  
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
