

kmer.compare = function(repA, repB, kmer.size = 8)
{
  if(is.na(repA))
  {
    print("repA not a character!")
    return(0)
  }
  if(is.na(repB))
  {
    print("repB not a character!")
    return(0)
  }
  repAkmers = extract_kmers(repA, kmer.size)
  repBkmers = extract_kmers(repB, kmer.size)
  common = ( length(which(repBkmers %in% repAkmers)) + length(which(repAkmers %in% repBkmers)) ) / 2
  if(is.na(common))
  {
    return(0)
  }
  d = 100 * common/ave(nchar(repA), nchar(repB))
  if(is.na(d))
  {
    return(0)
  }
  return(d)
}
