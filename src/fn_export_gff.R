
export.gff = function(temp.folder = execution.path, 
                       assemblyName = sequences$file.name[i])
{
  print("Export gff function")

  repeats = read.csv(file = paste(execution.path, "/all.repeats.from.", assemblyName, ".csv", sep = ""))
  
  
  gff_format = data.frame(seqid = vector(mode = "numeric", length = nrow(repeats)),
                          source = vector(mode = "numeric", length = nrow(repeats)),
                          type = vector(mode = "numeric", length = nrow(repeats)),
                          start = vector(mode = "numeric", length = nrow(repeats)),
                          end = vector(mode = "numeric", length = nrow(repeats)),
                          score = vector(mode = "numeric", length = nrow(repeats)),
                          strand = vector(mode = "numeric", length = nrow(repeats)),
                          phase = vector(mode = "numeric", length = nrow(repeats)),
                          attributes = vector(mode = "numeric", length = nrow(repeats)))
  
  
  gff_format$source = "TRASH"
  gff_format$type = "satellite_DNA"
  gff_format$start = repeats$start
  gff_format$end = repeats$end
  gff_format$score = "."
  gff_format$strand = repeats$strand
  gff_format$phase = "."
  
  
  for(i in 1 : nrow(gff_format))
  {
    print(i)
    gff_format$seqid[i] = strsplit(repeats$region.name[i], split = paste(assemblyName, "_", sep = ""))[[1]][2]
    gff_format$attributes[i] = paste("Class=", repeats$class[i], sep = "")
    
  }
  options(scipen=10)
  write.table(x = gff_format, file = paste(temp.folder, "/gff_", assemblyName, sep = ""), quote = FALSE, sep = "\t", eol = "\r", row.names = FALSE, col.names = FALSE)
  options(scipen=0)
  
}