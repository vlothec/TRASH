###### extract N bases with highest frequencies in the alignment

consensus_N = function(alignment, N)
{
  alignment.matrix = as.matrix.alignment(alignment)
  
  frequencies = vector(mode = "numeric", length = ncol(alignment.matrix))
  
  for(i in 1 : length(frequencies))
  {
    frequencies[i] = frequencies[i] + sum(1 * (alignment.matrix[,i] != "-"))
  }
  
  is.in.consensus = vector(mode = "logical", length = ncol(alignment.matrix))
  
  for(i in 1 : N)
  {
    is.in.consensus[order(frequencies, decreasing = TRUE)[i]] = TRUE
  }
  
  consensus = vector(mode = "logical", length = N)
  consensus.ID = 1
  for(i in 1 : ncol(alignment.matrix))
  {
    if(is.in.consensus[i])
    {
      consensus[consensus.ID] = c("g","c","t","a")[which.max(c(length(which(alignment.matrix[,i] == "g")),
                                                               length(which(alignment.matrix[,i] == "c")),
                                                               length(which(alignment.matrix[,i] == "t")),
                                                               length(which(alignment.matrix[,i] == "a"))))]
      consensus.ID = consensus.ID + 1
    }
  }
  
  return(toupper(paste(consensus, collapse = "")))
}





























