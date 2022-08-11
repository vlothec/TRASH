

extract_kmers = function(fasta.sequence, kmer.size)
{
  fasta.sequence = strsplit(fasta.sequence, split = "")[[1]]
  size = length(fasta.sequence)
  fasta.sequence = c(fasta.sequence, fasta.sequence)
  kmers = vector(mode = "character", length = size)
  for(i in 1 : size)
  {
    kmers[i] = paste(fasta.sequence[i : (i + kmer.size - 1)], collapse = "")
  }
  return(kmers)
}
