new.distance.N = function(plot = F, distances, N.max.div = 100, try.until = 12, smooth.percent = 2)
{
  N = NULL
  
  dist.hist = hist(distances[!is.na(distances)], breaks = max(distances[!is.na(distances)]), plot = FALSE) #smoothen by reducing breaks (need to do it consiously)
  
  #make the divisions
  
  
  N = dist.hist$breaks[which.max(dist.hist$counts) + 1]
  
  hist.smooth.range = ceiling(N*smooth.percent/100)
  
  new.values = vector(mode = "numeric", length = try.until)
  
  for(i in 1 : try.until)
  {
    new.values[i] = sum(dist.hist$counts[dist.hist$breaks %in% (ceiling((N-hist.smooth.range)/i) : ceiling((N+hist.smooth.range)/i)) ])
  }
  scores.with.multiplications = vector(mode = "numeric", length = try.until)
  for(i in 1 : try.until)
  {
    new.N = round(N/i)
    new.value.score = new.values[i]
    possible.multiplications = 1:i
    scores.with.multiplications[i] = sum(new.values[possible.multiplications])
    
  }
  if(plot)
  {
    #plot(scores.with.multiplications, ylim = c(0, max(scores.with.multiplications)))
    plot(dist.hist, xlim = c(0, 500))
  }
  
  options(scipen=999)
  scores.ratio = 100 * (1 - (scores.with.multiplications[1:(try.until-1)] / scores.with.multiplications[2:try.until]))
  
  if(length(which(scores.ratio > N.max.div)) > 0)
  {
    N = N / round(max(which(scores.ratio > N.max.div)) + 1)
  }
  
  return(N)
}









