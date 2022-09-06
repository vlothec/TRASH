
HOR.wrapper = function(threshold = 5, 
                       cutoff = 2, 
                       temp.folder = "", #
                       assemblyName = "", #sequences$file.name[i], 
                       chr.name = "",
                       mafft.bat.file = "", #paste(installation.path, "/src/mafft-linux64/mafft.bat", sep = ""),
                       hor.c.script.path = "",
                       class.name = "")
{
  print(paste("HOR function on ", chr.name, ", using repeats of class ", class.name, sep = ""))
  repeats = read.csv(file = paste(execution.path, "/all.repeats.from.", assemblyName, ".csv", sep = ""))
  classes = unique(repeats$class)
  
  print(paste("Available classes: ", classes, sep = ""))

  
  if(length(classes) == 0)
  {
    print(paste("no classes assigned to the repeats", sep = ""))
    return(1)
  }
  
  repeats$class[is.na(repeats$class)] = ""
    
  repeats = repeats[repeats$class == class.name,]

  
  if(nrow(repeats) == 0)
  {
    print(paste("no repeats of the class ", class.name, sep = ""))
    return(1)
  }
  
  chr.names = unique(repeats$region.name)
  
  repeats$region.name[is.na(repeats$region.name)] = ""
  
  repeats = repeats[repeats$region.name == paste(assemblyName, chr.name, sep = "_"),]
  
  if(nrow(repeats) == 0)
  {
    print(paste("no repeats of the name ", chr.name, " and class ", class.name, sep = ""))
    return(1)
  }
  
  work.dir = paste(temp.folder, "/HOR", sep = "")
  
  if(!dir.exists(work.dir))
  {
    dir.create(work.dir)
  }
  setwd(work.dir)
  
  
  repeats$strand[repeats$strand == "+"] = "1"
  repeats$strand[repeats$strand == "-"] = "2"
    
  file.to.align = paste("to_align_", assemblyName, "_", chr.name, "_", class.name, ".fasta", sep = "")
  
  write.fasta(sequences = as.list(repeats$seq), names = as.list(paste(1:length(repeats$strand), "_D", repeats$strand,  sep = "")), file.out = file.to.align, open = "w")
  
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
  
  output.alignment = paste("aligned_", assemblyName, "_", chr.name, "_", class.name, ".fasta", sep = "")
  
  system(paste(mafft.bat.file, " --quiet --inputorder --kimura 1 --retree 1 ", file.to.align, " > ", output.alignment, sep = ""), intern = TRUE)
  
  if(!file.exists(output.alignment))
  {
    print("cannot find alignment file for HORs, stopping HOR calculation")
    return(1)
  }
  
  if(file.size(output.alignment) == 0)
  {
    print("output alignment file for HORs empty, stopping HOR calculation")
    return(1)
  }
  
  print("Running HOR software")
  
  system(paste(hor.c.script.path, " ", work.dir, "/ ", output.alignment, " ", 1, " ", threshold, " ", cutoff, " ", 1, sep = ""), intern = FALSE)
  
  hor.output.file = paste(work.dir, "/", "HORs_method_1_", output.alignment, "_t_", threshold, "_c_", cutoff, ".csv", sep = "")
  
  if(!file.exists(hor.output.file))
  {
    print("cannot find output file for HORs, stopping HOR calculation")
    return(1)
  }
  
  if(file.size(hor.output.file) == 0)
  {
    print("output HOR file empty, stopping HOR calculation")
    return(1)
  }
  
  hors = read.csv(file = hor.output.file, header = TRUE)
  
  if(nrow(hors) > 1)
  {
    hors = as.data.frame(hors)
    
    hors$start.A.bp = repeats$start[hors$start_A]
    hors$start.B.bp = repeats$start[hors$start_B]
    hors$end.A.bp = repeats$end[hors$end_A]
    hors$end.B.bp = repeats$end[hors$end_B] 
    hors$chrA = repeats$region.name[hors$start_A]  
    hors$chrB = repeats$region.name[hors$start_B]
    hors$block.size.in.units = hors$end_A - hors$start_A + 1
    hors$block.A.size.bp = hors$end.A.bp - hors$start.A.bp + 1
    hors$block.B.size.bp = hors$end.B.bp - hors$start.B.bp + 1
    for(p in nrow(hors) : 1)
    {
      block1size = abs(hors$end.A.bp[p] - hors$start.A.bp[p])
      block2size = abs(hors$end.B.bp[p] - hors$start.B.bp[p])
      if(abs(block2size - block1size) > 1000)
      {
        hors = hors[-p,]
        p = p - 1
      }
    }
    
    write.csv(x = hors, file = paste(temp.folder, "/HORs_", assemblyName, "_", chr.name, "_", class.name, ".csv", sep = ""))
    
    png(filename = paste(work.dir, "/HOR_plot_", assemblyName, "_", chr.name, "_", class.name, ".png", sep = ""), width = 5000, height = 5000, pointsize = 25)
    plot(xlab = paste(assemblyName, "_", chr.name, "_", class.name, sep = ""), 
         ylab = paste(assemblyName, "_", chr.name, "_", class.name, sep = ""), 
         x = hors$start.A.bp, pch = 19, cex = 0.1, 
         y = hors$start.B.bp, 
         xlim = c(min(hors$start.A.bp), max(hors$start.A.bp)), 
         ylim = c(min(hors$start.B.bp), max(hors$start.B.bp)))
    dev.off()
  }
    
  print("HOR done")
  return(0)
}
