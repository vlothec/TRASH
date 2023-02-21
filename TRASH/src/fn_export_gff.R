edit.seq.id = function(region.name.string = "", assembly.name = "")
{
  return(strsplit(region.name.string, split = paste(assembly.name, "_", sep = ""))[[1]][2])
}
edit.attributes = function(class.string = "")
{
  return(paste("Class=", class.string, sep = ""))
}

export.gff = function(temp.folder = "", 
                       assemblyName = "")
{
  print("Export gff function")

  repeats = read.csv(file = paste(execution.path, "/all.repeats.from.", assemblyName, ".csv", sep = ""))
  
  print("read repeats")
  
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
  
  gff_format$seqid = unlist(lapply(X = repeats$region.name, edit.seq.id, assemblyName))
  gff_format$attributes = unlist(lapply(X = repeats$class, edit.attributes))
  
  
  print("wrtiting repeats")
  options(scipen=10)
  write.table(x = gff_format, file = paste(temp.folder, "/gff_", assemblyName, sep = ""), quote = FALSE, sep = "\t", eol = "\r", row.names = FALSE, col.names = FALSE)
  options(scipen=0)
  
}