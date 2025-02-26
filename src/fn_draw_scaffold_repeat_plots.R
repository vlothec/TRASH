

draw.scaffold.repeat.plots = function(temp.folder = "", assemblyName = "", fastaDirectory = "", single.pngs = FALSE, y.max = 600)
{
  fasta = read.fasta(fastaDirectory, as.string = TRUE)
  if(temp.folder == "")
  {
    stop("no temp folder")
  }
  setwd(temp.folder)
  
  if(!dir.exists(paste(assemblyName, "_out", sep = "")))
  {
    dir.create(paste(assemblyName, "_out", sep = ""))
  }
  
  if(!dir.exists(paste(temp.folder, "/", paste(assemblyName, "_out", sep = ""), "/plots", sep = "")))
  {
    dir.create(paste(temp.folder, "/", paste(assemblyName, "_out", sep = ""), "/plots", sep = ""))
  }
  
  
  
  files = list.files(path = paste(temp.folder, "/", paste(assemblyName, "_out", sep = ""), "/temp2", sep = ""), pattern = "Repeats*")
  
  rbPal = colorRampPalette(c('red','blue'))
  
  if(single.pngs == FALSE)
  {
    #save as a single pdf
    i = 1
    pdf(file = paste(temp.folder, "/plots/", assemblyName, "_linearplot.pdf", sep = ""), width = 50, height = 10, pointsize = 30)
    while(i < length(files))
    {
      k = i
      scaffold = str_split(files[i], pattern = "_")[[1]]
      scaffold = paste(scaffold[2 : (length(scaffold) - 2)], collapse = "_")
      fragment.length = as.numeric(nchar(fasta[which(names(fasta) == scaffold)]))
      if(length(fragment.length) != 0)
      {
        plot(0,0, type="p", xlab="repeat start position", ylab="repeat size", xlim=c(0, fragment.length), ylim=c(0, y.max), main = scaffold, cex = 0)
        
        while(k <= length(files) && paste((str_split(files[k], pattern = "_")[[1]])[2 : (length((str_split(files[k], pattern = "_")[[1]])) - 2)], collapse = "_") == scaffold)
        {
          repeats = read.csv(file = paste(temp.folder, "/", paste(assemblyName, "_out", sep = ""), "/temp2/", files[k], sep = ""), header = TRUE)
          if(nrow(repeats) > 0)
          {
            repeats$start = repeats$start + as.numeric(str_split(files[k], pattern = "_")[[1]][length(str_split(files[k], pattern = "_")[[1]]) - 1])
            repeats$col = rbPal(y.max)[repeats$width]
            points(repeats$start,repeats$width, col = repeats$col, pch = ".", cex = 7)
          }
          remove(repeats)
          k = k + 1
        }
      }
      
      i = k
    }
    dev.off()
  }
  if(single.pngs == TRUE)
  {
    i = 1
    while(i < length(files))
    {
      k = i
      scaffold = str_split(files[i], pattern = "_")[[1]]
      scaffold = paste(scaffold[2 : (length(scaffold) - 2)], collapse = "_")
      fragment.length = as.numeric(nchar(fasta[which(names(fasta) == scaffold)]))
      
      if(length(fragment.length) > 0)
      {
        png(filename = paste(temp.folder, "/plots/", assemblyName, "_", i, ".png", sep = ""), width = 4000, height = 1000, pointsize = 45)
        plot(0,0, type="p", xlab="repeat start position", ylab="repeat size", xlim=c(0, fragment.length), ylim=c(0, y.max), main = scaffold, cex = 0)
        
        while(k <= length(files) && paste((str_split(files[k], pattern = "_")[[1]])[2 : (length((str_split(files[k], pattern = "_")[[1]])) - 2)], collapse = "_") == scaffold)
        {
          repeats = read.csv(file = paste(temp.folder, "/", paste(assemblyName, "_out", sep = ""), "/temp2/", files[k], sep = ""), header = TRUE)
          if(nrow(repeats) > 0)
          {
            repeats$start = repeats$start + as.numeric(str_split(files[k], pattern = "_")[[1]][length(str_split(files[k], pattern = "_")[[1]]) - 1])
            repeats$col = rbPal(y.max)[repeats$width]
            points(repeats$start,repeats$width, col = repeats$col, pch = ".", cex = 7)
          }
          remove(repeats)
          k = k + 1
        }
        dev.off()
      }
      i = k
    }
  }
  return(0)
}
