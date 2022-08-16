Repeat.Identifier = function(DNA.sequence = "", assemblyName = "", fasta.name = "", 
                             kmer = 12, window = 1000, threshold = 10, mask.small.regions = 1500, mask.small.repeats = 4,
                             max.repeat.size = 500,
                             tests = 5, temp.folder = "",
                             sequence.template, mafft.bat.file = "", LIMIT.REPEATS.TO.ALIGN = 78000)
{
  
  #TODO check for all the parameters given into the function that are correct here
  
  print("Repeat Identification new")
  if(temp.folder == "")
  {
    stop("no temp folder")
  }
  print(temp.folder)
  print(assemblyName)
  setwd(temp.folder)
  
  if(!dir.exists(paste(assemblyName, "_out", sep = "")))
  {
    dir.create(paste(assemblyName, "_out", sep = ""))
  }
  
  write(paste("starting the main function", assemblyName, fasta.name, sep = " "), file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
  
  regions.data.frame = data.frame(index = integer(), fasta.name = character(), start = integer(), end = integer(), ave.score = double(), most.freq.value.N = integer())#####
  
  
  seq.length = nchar(DNA.sequence)
  
  #calculate kmer scores for each window
  {
    scores = vector(mode = "numeric", length = seq.length %/% window + 1)
    for(ii in 1 : (length(scores) - 1))  # for every window of size window but not the last one calculate kmer score
    {
      if(ii%%100 == 0){
        write(paste("Checking windows of sequence ", fasta.name, ", window ", ii, "/", length(scores), sep = ""), file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
      }
      window.start = (ii - 1) * window
      check.window = str_sub(DNA.sequence, window.start, (window.start + window - 1))
      for(iii in 1 : (window - kmer))
      {
        check.string = str_sub(string = check.window, start = iii, end = (iii + kmer - 1))
        count = 1 * ((str_count(string = check.string, pattern = "N") == 0) & str_count(string = check.window, pattern = check.string) > 1)
        scores[ii] = scores[ii] + count
      }
    }
    ii = length(scores) # for the last window
    write(paste("sequence ", fasta.name, ", window ", ii, "/", length(scores), sep = ""), file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
    window.start = (ii - 1) * window
    check.window = str_sub(DNA.sequence, window.start)
    if(nchar(check.window) > kmer)
    {
      for(iii in 1 : (nchar(check.window) - kmer))
      {
        check.string = str_sub(string = check.window, start = iii, end = (iii + kmer - 1))
        count = 1 * ((str_count(string = check.string, pattern = "N") == 0) & str_count(string = check.window, pattern = check.string) > 1)
        scores[ii] = scores[ii] + count
      }
    }
    
    scores[1 : (length(scores) - 1)] = scores[1 : (length(scores) - 1)] / window * 100 #change score into a percentage score
    scores[length(scores)] = scores[length(scores)] / nchar(check.window)
  }
  write(paste("sequence ", fasta.name, " finished scoring windows", sep = ""), file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
  
  #write(paste("sequence ", fasta.name, " scores: ", scores, sep = ""), file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
  
  
  #get continuous regions
  {
    smoothing.window.scores.across = 2
    max.sd.window.scores = 2
    index = 0
    ii = 1
    while(ii <= length(scores)) #apply threshold and save repetitive windows from this fasta
    {
      start.index = ii
      if(scores[ii] >= threshold)
      {
        index = index + 1
        name = fasta.name
        start = (ii - 1) * window
        end = start - 1
        ave.score = 0
        s = 0
        while(scores[ii] >= threshold)
        {
          s = s + 1
          end = end + window
          ave.score = scores[ii] + ave.score 
          ii = ii + 1
        }
        ave.score = ave.score / s
        most.freq.value.N = 0
        temp = data.frame(index, name, start, end, ave.score, most.freq.value.N)
        regions.data.frame = rbind(regions.data.frame, temp)
      }
      
      #check the windows for continuity and divide if necessary
      {
        if(FALSE)
        {
          #TODO NEW THING DONT DO
          if((ii - start) >= (smoothing.window.scores.across * 2))
          {
            averages.vector = vector(mode = "numeric", length = (ii - start.index + 1))
            sd.vector = vector(mode = "numeric", length = (ii - start.index + 1))
            ID.averages = 1
            for(j in start.index : ii)
            {
              indexes.to.check = (j - smoothing.window.scores.across) : (j + smoothing.window.scores.across)
              indexes.to.check = indexes.to.check[indexes.to.check > 0 & indexes.to.check < ii]
              averages.vector[ID.averages] = mean(scores[indexes.to.check])
              sd.vector[ID.averages] = sd(scores[indexes.to.check])
              ID.averages = ID.averages + 1
            }
            k = 1
            ranges.sd.high = which(sd.vector > max.sd.window.scores)
            while(k < length(ranges.sd.high))
            {
              if(ranges.sd.high[k+1] - ranges.sd.high[k] == 1)
                k = k + 1
            }
          }
        }
      }
      
      
      ii = ii + 1
    }
  }
  
  write(paste("sequence ", fasta.name, "finished calculating regions", sep = ""), file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
  
  #remove small regions
  regions.data.frame = regions.data.frame[(regions.data.frame$end - regions.data.frame$start) > mask.small.regions,]
  
  
  
  
  
  
  #stop if nothing identified
  if(nrow(regions.data.frame) < 1){
    write(paste("sequence ", fasta.name, "no repeats found", sep = ""), file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
    
    return(0)
  }
  
  write(paste("Step 1 (identify windows) done", sep = ""), file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
  
  #calculate distances between kmers from regions and find most common value N
  {
    for(i in 1 : nrow(regions.data.frame))
    {
      seqB = str_sub(DNA.sequence, regions.data.frame$start[i], regions.data.frame$end[i])
      write(paste("calculating distances on region ", i, "/", nrow(regions.data.frame), " region size: ", nchar(seqB), sep = ""), file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
      
      distance = vector(mode = "numeric", length =  nchar(seqB))
      for(ii in 1 : (nchar(seqB) - kmer + 1))
      {
        kmer.pattern = str_sub(seqB, ii, (ii + kmer - 1))
        window.string = str_sub(seqB, (ii + 1 + mask.small.repeats), (ii + max.repeat.size))
        distance[ii] = str_locate(string = window.string, pattern = kmer.pattern)[[1]] + mask.small.repeats
      }
      distance = distance[!is.na(distance)]
      distance = distance[distance > mask.small.repeats]
      if(length(distance) > 0)
      {
        hist.values = hist(distance, breaks = max(distance), plot = FALSE)
        regions.data.frame$most.freq.value.N[i] = hist.values$breaks[which.max(hist.values$counts) + 1]
      }
    }
    regions.data.frame = regions.data.frame[regions.data.frame$most.freq.value.N >= mask.small.repeats,]
  }
  write(paste("sequence ", fasta.name, ": finished calculating main distances", sep = ""), file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
  
  
  
  #stop if nothing identified
  if(nrow(regions.data.frame) < 1){
    
    write(paste("sequence ", fasta.name, "no regions identified", sep = ""), file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
    
    return(0)
  }
  if(!dir.exists(paste(temp.folder, "/", paste(assemblyName, "_out", sep = ""), "/temp1", sep ="")))
  {
    dir.create(paste(temp.folder, "/", paste(assemblyName, "_out", sep = ""), "/temp1", sep =""))
  }
  write.csv(regions.data.frame, file = paste(paste(assemblyName, "_out", sep = ""), "/temp1/Pre_Primary_regions_", assemblyName, "_", fasta.name, ".csv", sep = "")) #can silence
  
  #Test random N-size mers, handle overlaps, choose the best one, extract repeats for mafft alignment and align to extract primary consensus
  {
    regions.data.frame$consensus.primary = ""
    regions.data.frame$consensus.count = ""
    
    #for each region find primary consensus
    i = 1
    while(i <= nrow(regions.data.frame))
    {
      if(regions.data.frame$consensus.primary[i] == "")
      {
        write(paste("testing random N samples from region ", i, "/", nrow(regions.data.frame), " window size: ", regions.data.frame$end[i] - regions.data.frame$start[i], sep = ""), 
              file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
        seqC = str_sub(DNA.sequence, regions.data.frame$start[i], regions.data.frame$end[i])
        N = regions.data.frame$most.freq.value.N[i]
        if(N == 0)
        {
          N = 1
        }
        random.sequences.start = sample(1 : (nchar(seqC) - N), tests)
        random.sequences = str_sub(seqC, random.sequences.start, (random.sequences.start + N - 1))
        random.sequence.scores = vector(mode = "numeric", length = tests)
        random.sequence.counts = vector(mode = "numeric", length = tests)
        
        maxMis = N %/% 3
        if(maxMis > 100)
        {
          maxMis = 100
        }
        #test a number of random substrings to find best representation, remove overlapping 
        for(ii in 1 : tests)
        {
          match = matchPattern(pattern = random.sequences[ii], subject = seqC, max.mismatch = maxMis, with.indels = TRUE)
          temp = as.data.frame(match)
          if(nrow(temp) > 0)
          {
            temp$start  = start(match)
            temp$end = end(match)
            temp$strand = "+"
            iii = nrow(temp)
            random.sequence.counts[ii] = 0
            if(iii > 0)
            {
              random.sequence.counts[ii] = nrow(temp) 
            }
            while(iii > 1)
            {
              if((temp$start[iii] + 2 * maxMis) <  temp$end[iii - 1])
              {
                temp = temp[-iii,]
              }
              iii = iii - 1
            }
            if(nrow(temp) > 1)
            {
              write.csv(temp, file = paste(temp.folder, "/", paste(assemblyName, "_out", sep = ""), "/temp1/", fasta.name, "_", ii, ".csv", sep = ""), row.names = FALSE)
              random.sequence.scores[ii] = sum(temp$end - temp$start)
            }
          }
        }
        
        #mafft the best representation and extract primary consensus
        for(L in 1 : length(random.sequence.scores))
        {
          if(!is.na(random.sequence.scores[L])) #if any of the scores is good, continue and get out of the loop, otherwise there will be nothing assigned
          {
            if(random.sequence.scores[which.max(random.sequence.scores)] > 0)
            {
              winner = read.csv(paste(temp.folder, "/", paste(assemblyName, "_out", sep = ""), "/temp1/", fasta.name, "_", which.max(random.sequence.scores), ".csv", sep = ""))
              
              if(nrow(winner) > 0)
              {
                if((nrow(winner) * N) > LIMIT.REPEATS.TO.ALIGN)
                {
                  sample.IDs = sample(x = 1:nrow(winner), size = round(LIMIT.REPEATS.TO.ALIGN/N, 0))
                } else
                {
                  sample.IDs = 1:nrow(winner)
                }
                if(ncol(winner) == 5)
                {
                  write.fasta(file.out = paste(paste(assemblyName, "_out", sep = ""), "/temp1/", fasta.name, "_pre.extract.fasta", sep = ""), 
                              names = paste(winner$start, winner$end, sep = "_")[sample.IDs], 
                              sequences = str_split(winner[,4], pattern = "")[sample.IDs],
                              as.string = FALSE)
                } else if(ncol(winner) == 4){
                  write.fasta(file.out = paste(paste(assemblyName, "_out", sep = ""), "/temp1/", fasta.name, "_pre.extract.fasta", sep = ""), 
                              names = paste(winner$start, winner$end, sep = "_")[sample.IDs], 
                              sequences = str_split(winner$x, pattern = "")[sample.IDs], 
                              as.string = FALSE)
                }
                input = paste(paste(temp.folder, "/", assemblyName, "_out", sep = ""), "/temp1/", fasta.name, "_pre.extract.fasta", sep = "")
                output = paste(paste(temp.folder, "/", assemblyName, "_out", sep = ""), "/temp1/", fasta.name, "_pre.extract_aligned.fasta", sep = "")
                
                system(paste(mafft.bat.file, " --quiet --retree 2 --inputorder ", input, " > ", output, sep = ""), intern = TRUE)
                sample.IDs = NULL
                
                alignment = read.alignment(output, format = "FASTA", forceToLower = FALSE)
                N = regions.data.frame$most.freq.value.N[i]
                if(nchar(alignment$seq[[1]]) < N)
                {
                  consensus = consensus_N(alignment, nchar(alignment$seq[[1]]))
                }
                else
                {
                  consensus = consensus_N(alignment, N)
                }
                
                
                regions.data.frame$consensus.primary[i] = consensus
                regions.data.frame$consensus.count[i] = nrow(winner)
                remove(alignment, input, output)
                
              } else
              {
                regions.data.frame$consensus.primary[i] = "none_identified"
                regions.data.frame$consensus.count[i] = 0
              }
            } else
            {
              regions.data.frame$consensus.primary[i] = "none_identified"
              regions.data.frame$consensus.count[i] = 0
            }
            break()
          }
        }
        if(regions.data.frame$consensus.primary[i] == "") #in case all scores were rubbish
        {
          regions.data.frame$consensus.primary[i] = "none_identified"
          regions.data.frame$consensus.count[i] = 0
        }
        
        #if there was some coverage, but not all, separate the identified region and make a new one to identify the next class of repeats
        if(regions.data.frame$consensus.primary[i] != "none_identified")
        {
          if(regions.data.frame$consensus.primary[i] != "")
          {
            
            new.regions.starts = NULL
            new.regions.ends = NULL
            
            
            max.score.for.region = nchar(seqC)
            identified.score = random.sequence.scores[which.max(random.sequence.scores)]
            if(is.numeric(identified.score) & is.numeric(max.score.for.region))
            {
              leftover.space = max.score.for.region - identified.score
              if(identified.score > mask.small.regions)#see if total identified space is more than minimum region size
              {
                if(leftover.space > mask.small.regions) #see if total leftover space is more than minimum region size
                {
                  write(paste("find empty space ", i, " with space left: ", leftover.space, sep = " "), file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
                  write.csv(regions.data.frame, file = paste(temp.folder, "/", assemblyName, "_out/temp1/Temp_regions_", paste(assemblyName, "_out", sep = ""), "_", fasta.name, ".csv", sep = ""))
                  
                  #find consecutive windows that would be unoccupied by the repeats and are bigger than minimum region size
                  
                  start.repetitive.region = regions.data.frame$start[i]
                  end.repetitive.region = regions.data.frame$end[i]
                  
                  start.new.region = start.repetitive.region
                  
                  if(nrow(winner) > 1)
                  {
                    for(s in 1 : nrow(winner))
                    {
                      end.new.region = winner$start[s] + start.repetitive.region
                      
                      if((end.new.region - start.new.region) > mask.small.regions) #make a new region, also find a new N value for it
                      {
                        new.distance.N = 0
                        new.index = nrow(regions.data.frame) + 1
                        
                        seqB = str_sub(DNA.sequence, start.new.region, end.new.region)
                        write(paste("calculating distances on region ", new.index, "/", (nrow(regions.data.frame) + 1), " region size: ", nchar(seqB), sep = ""), file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
                        
                        distance = vector(mode = "numeric", length =  nchar(seqB))
                        for(ii in 1 : (nchar(seqB) - kmer + 1))
                        {
                          kmer.pattern = str_sub(seqB, ii, (ii + kmer - 1))
                          window.string = str_sub(seqB, (ii + 1 + mask.small.repeats), (ii + max.repeat.size))
                          distance[ii] = str_locate(string = window.string, pattern = kmer.pattern)[[1]] + mask.small.repeats
                        }
                        distance = distance[!is.na(distance)]
                        distance = distance[distance > mask.small.repeats]
                        if(length(distance) > 0)
                        {
                          hist.values = hist(distance, breaks = max(distance), plot = FALSE)
                          new.distance.N = hist.values$breaks[which.max(hist.values$counts) + 1]
                        }
                        
                        if(!is.na(new.distance.N) & !is.null(new.distance.N) & new.distance.N > 0)
                        {
                          write(paste("new row: ", new.index, 
                                      fasta.name, 
                                      start.new.region,
                                      end.new.region, 
                                      regions.data.frame$ave.score[i], 
                                      new.distance.N, "", 0, sep = ", "),
                                file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
                          
                          regions.data.frame[(nrow(regions.data.frame) + 1),] = c((max(regions.data.frame$index) + 1), 
                                                                                  fasta.name, 
                                                                                  start.new.region,
                                                                                  end.new.region, 
                                                                                  regions.data.frame$ave.score[i], 
                                                                                  new.distance.N, "", 0)
                          
                          new.regions.starts = c(new.regions.starts, as.numeric(start.new.region))
                          new.regions.ends = c(new.regions.ends, as.numeric(end.new.region))
                          
                        }
                        
                        start.new.region = winner$end[s] + start.repetitive.region
                      }
                      
                      start.new.region = winner$end[s] + start.repetitive.region
                    }
                    remove(winner)
                  }
                  if(start.new.region != start.repetitive.region)
                  {
                    end.new.region = end.repetitive.region
                    
                    #for the last time after the final repeat
                    write(paste("find empty space last time", i, " with space left: ", leftover.space, sep = " "), file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
                    
                    if((end.new.region - start.new.region) > mask.small.regions) #make a new region, also find a new N value for it
                    {
                      new.distance.N = 0
                      new.index = nrow(regions.data.frame) + 1
                      
                      seqB = str_sub(DNA.sequence, start.new.region, end.new.region)
                      write(paste("calculating distances on region ", new.index, "/", (nrow(regions.data.frame) + 1), " region size: ", nchar(seqB), sep = ""), file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
                      
                      distance = vector(mode = "numeric", length =  nchar(seqB))
                      for(ii in 1 : (nchar(seqB) - kmer + 1))
                      {
                        kmer.pattern = str_sub(seqB, ii, (ii + kmer - 1))
                        window.string = str_sub(seqB, (ii + 1 + mask.small.repeats), (ii + max.repeat.size))
                        distance[ii] = str_locate(string = window.string, pattern = kmer.pattern)[[1]] + mask.small.repeats
                      }
                      distance = distance[!is.na(distance)]
                      distance = distance[distance > mask.small.repeats]
                      if(length(distance) > 0)
                      {
                        hist.values = hist(distance, breaks = max(distance), plot = FALSE)
                        new.distance.N = hist.values$breaks[which.max(hist.values$counts) + 1]
                      }
                      
                      if(!is.na(new.distance.N) & !is.null(new.distance.N) & new.distance.N > 0)
                      {
                        write(paste("new row: ", new.index, 
                                    fasta.name, 
                                    start.new.region,
                                    end.new.region, 
                                    regions.data.frame$ave.score[i], 
                                    new.distance.N, "", 0, sep = ", "),
                              file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
                        regions.data.frame[(nrow(regions.data.frame) + 1),] = c((max(regions.data.frame$index) + 1), 
                                                                                fasta.name, 
                                                                                start.new.region,
                                                                                end.new.region, 
                                                                                regions.data.frame$ave.score[i], 
                                                                                new.distance.N, "", 0)
                        
                        new.regions.starts = c(new.regions.starts, as.numeric(start.new.region))
                        new.regions.ends = c(new.regions.ends, as.numeric(end.new.region))
                      }
                    }
                  }
                }
              }
            }
          }
        }
        regions.data.frame$start = as.numeric(regions.data.frame$start)
        regions.data.frame$end = as.numeric(regions.data.frame$end)
        regions.data.frame$index = as.numeric(regions.data.frame$index)
        regions.data.frame$ave.score = as.numeric(regions.data.frame$ave.score)
        regions.data.frame$most.freq.value.N = as.numeric(regions.data.frame$most.freq.value.N)
        regions.data.frame$consensus.count = as.numeric(regions.data.frame$consensus.count)
        
        
        
        if(length(new.regions.starts) > 0)
        {
          old.start = start.repetitive.region
          for(d in 1 : length(new.regions.starts))
          {
            old.end = new.regions.starts[d]
            if(old.end > old.start)
            {
              write(paste("new old region split d: ", d, sep = ""), file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
              
              regions.data.frame[(nrow(regions.data.frame) + 1),] = c(as.numeric(max(regions.data.frame$index) + 1), 
                                                                      fasta.name, 
                                                                      as.numeric(old.start),
                                                                      as.numeric(old.end), 
                                                                      as.numeric(regions.data.frame$ave.score[i]), 
                                                                      as.numeric(regions.data.frame$most.freq.value.N[i]), 
                                                                      regions.data.frame$consensus.primary[i], 
                                                                      0)
              regions.data.frame$name[i] = "REMOVE_ROW"
            }
            old.start = new.regions.ends[d]
          }
          old.end = end.repetitive.region #for the last one
          if(old.end > old.start)
          {
            write(paste("new old region split d: ", d, sep = ""), file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
            
            regions.data.frame[(nrow(regions.data.frame) + 1),] = c(as.numeric(max(regions.data.frame$index) + 1), 
                                                                    fasta.name, 
                                                                    as.numeric(old.start),
                                                                    as.numeric(old.end), 
                                                                    as.numeric(regions.data.frame$ave.score[i]), 
                                                                    as.numeric(regions.data.frame$most.freq.value.N[i]), 
                                                                    regions.data.frame$consensus.primary[i], 
                                                                    0)
            regions.data.frame$name[i] = "REMOVE_ROW"
          }
        }
        regions.data.frame$start = as.numeric(regions.data.frame$start)
        regions.data.frame$end = as.numeric(regions.data.frame$end)
        regions.data.frame$index = as.numeric(regions.data.frame$index)
        regions.data.frame$ave.score = as.numeric(regions.data.frame$ave.score)
        regions.data.frame$most.freq.value.N = as.numeric(regions.data.frame$most.freq.value.N)
        regions.data.frame$consensus.count = as.numeric(regions.data.frame$consensus.count)
        
        if(i > 10000)
        {
          print("TRASH warning: number of identified regions surpassed 10,000")
        }
        write(paste("next loop i: ", i, " nrow data = ", nrow(regions.data.frame), sep = ""), file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
        
        write.csv(regions.data.frame, file = paste(temp.folder, "/", assemblyName, "_out/temp1/Temp_regions_", paste(assemblyName, "_out", sep = ""), "_", fasta.name, ".csv", sep = ""))
      }
      i = i + 1
    }
    write(paste("remove regions id: ", which(regions.data.frame$name == "REMOVE_ROW"), sep = ""), file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
    
    if(length(which(regions.data.frame$name == "REMOVE_ROW")) > 0)
    {
      regions.data.frame = regions.data.frame[(regions.data.frame$name != "REMOVE_ROW"),]
    }
    
    
    write(paste("out of the primary consensus loop", sep = ""), file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
    
    
    
    write(paste("finished test random Nmers of ", assemblyName, ": ", fasta.name, sep = ""), file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
    if(nrow(regions.data.frame) < 1)
    {
      write(paste("sequence ", fasta.name, " no regions after filtering", sep = ""), file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
      
      return(0)
    }
    write.csv(regions.data.frame, file = paste(temp.folder, "/", assemblyName, "_out/temp1/Pre_shift_regions_", paste(assemblyName, "_out", sep = ""), "_", fasta.name, ".csv", sep = ""))
    
    regions.data.frame$fasta.name = as.character(fasta.name)
  }
  
  #hash each primary consensus
  for(i in 1 : nrow(regions.data.frame))
  {
    if(!is.na(regions.data.frame$consensus.primary[i]))
    {
      if(regions.data.frame$consensus.primary[i] != "none_identified")
      {
        regions.data.frame$consensus.primary[i] = Hash_And_Reverse(regions.data.frame$consensus.primary[i], 3)
      }
    }
  }
  
  #shift the primary consensus sequences to match a given sequence, only those regions with N within a given range
  regions.data.frame$class = ""
  if(!is.na(sequence.template)[1])
  {
    for(i in 1 : nrow(regions.data.frame))
    {
      if(!is.na(regions.data.frame$consensus.primary[i]))
      {
        if(regions.data.frame$consensus.primary[i] != "none_identified")
        {
          highest = 0
          check.seq = regions.data.frame$consensus.primary[i]
          if((nchar(check.seq) > (min(sequence.template$length)-10)) & (nchar(check.seq) < (max(sequence.template$length)+10)))
          {
            highest = 0
            scores = 0
            for(ii in 1 : nrow(sequence.template))
            {
              if((nchar(check.seq) > (sequence.template$length[ii])*0.9) & (nchar(check.seq) < (sequence.template$length[ii])*1.1))
              {
                a = kmer.compare(sequence.template$seq[ii], check.seq)
                if(a > scores)
                {
                  scores = a
                  highest = ii
                }
                a = kmer.compare(sequence.template$seq[ii], revCompString(check.seq))
                if(a > scores)
                {
                  scores = a
                  highest = ii
                }
              }
            }
            if(highest != 0)
            {
              seq.len = nchar(check.seq)
              #create a df of all possible shifts and a column for their scores
              shifts = data.frame(seq = vector(mode = "character", length = (seq.len * 2)), score = vector(mode = "numeric", length = (seq.len * 2)))
              #split the sequence to make it easier to handle
              seq.to.split.plus = strsplit(check.seq, split = "")[[1]]
              seq.to.split.revcomp = strsplit(revCompString(check.seq), split = "")[[1]]
              #handle problematic first and last cases
              shifts$seq[1] = paste(seq.to.split.plus, collapse = "")
              shifts$seq[(seq.len)] = paste(c(seq.to.split.plus[seq.len], seq.to.split.plus[1 : (seq.len - 1)]), collapse = "")
              shifts$seq[seq.len  + 1] = paste(seq.to.split.revcomp, collapse = "")
              shifts$seq[(seq.len * 2)] = paste(c(seq.to.split.revcomp[seq.len], seq.to.split.revcomp[1 : (seq.len - 1)]), collapse = "")
              #do forward strand
              for(ii in 2 : (seq.len - 1))
              {
                shifts$seq[ii] = paste(c(seq.to.split.plus[ii : seq.len], seq.to.split.plus[1 : (ii - 1)]), collapse = "")
              }
              #do reverse strand
              for(ii in 2 : (seq.len - 1))
              {
                shifts$seq[ii + seq.len] = paste(c(seq.to.split.revcomp[ii : seq.len], seq.to.split.revcomp[1 : (ii - 1)]), collapse = "")
              }
              #align sequence.template to the shifts
              for(ii in 1 : nrow(shifts))
              {
                if(ii%%100 == 0)
                {
                  write(paste("Aligning shfit ", ii, " out of ", nrow(shifts), " from region ", i, paste = ""), file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
                }
                #TODO choose an alignment method that gives a precise score and align and save the score in the table
                shifts$score[ii] =  pairwiseAlignment(pattern = shifts$seq[ii], subject = sequence.template$seq[highest], type = "global", scoreOnly = TRUE)
              }
              #return the highest score shift as a primary.shifted[i]
              regions.data.frame$consensus.primary[i] = shifts$seq[which.max(shifts$score)]
              regions.data.frame$class[i] = sequence.template$name[highest]
            }
          }
        }
      }
    }
  }
  
  
  
  
  #Refine each consensus by mapping again, aligning with mafft and extracting consensus 
  {
    if(temp.folder == "")
    {
      stop("no temp folder")
    }
    setwd(temp.folder)
    
    if(!dir.exists(paste(temp.folder, "/", assemblyName, "_out", sep = "")))
    {
      dir.create(paste(temp.folder, "/", assemblyName, "_out", sep = ""))
    }
    if(!dir.exists(paste(paste(temp.folder, "/", assemblyName, "_out", sep = ""), "/temp2", sep ="")))
    {
      dir.create(paste(paste(temp.folder, "/", assemblyName, "_out", sep = ""), "/temp2", sep =""))
    }
    if(!dir.exists(paste(paste(temp.folder, "/", assemblyName, "_out", sep = ""), "/temp2/inputs", sep = "")))
    {
      dir.create(paste(paste(temp.folder, "/", assemblyName, "_out", sep = ""), "/temp2/inputs", sep = ""))
    }
    if(!dir.exists(paste(paste(temp.folder, "/", assemblyName, "_out", sep = ""), "/temp2/outputs", sep = "")))
    {
      dir.create(paste(paste(temp.folder, "/", assemblyName, "_out", sep = ""), "/temp2/outputs", sep = ""))
    }
    
    consensus.secondary = vector(mode = "character", length = nrow(regions.data.frame))
    repeats.identified = vector(mode = "numeric", length = nrow(regions.data.frame))
    
    
    for(i in 1 : nrow(regions.data.frame))
    {
      if(regions.data.frame$consensus.primary[i] != "none_identified")
      {
        if(nchar(regions.data.frame$consensus.primary[i]) > mask.small.repeats)
        {
          write(x =  paste("on ", regions.data.frame$fasta.name[i], ", generating consensus for region ", i, "/", nrow(regions.data.frame), " window size: ", regions.data.frame$end[i] - regions.data.frame$start[i], sep = ""), 
                file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
          seqC = str_sub(DNA.sequence, regions.data.frame$start[i], regions.data.frame$end[i])
          N = regions.data.frame$most.freq.value.N[i]
          maxMis = (N %/% 3) - 2
          if(maxMis > 100)
          {
            maxMis = 100
          } else if(maxMis < 1)
          {
            maxMis = 1
          }
          matchPlus = matchPattern(pattern = regions.data.frame$consensus.primary[i], subject = seqC, max.mismatch = maxMis, with.indels = TRUE)
          tempP = as.data.frame(matchPlus)
          tempP$start  = start(matchPlus)
          tempP$end = end(matchPlus)
          if(nrow(tempP) > 0)
          {
            tempP$strand = "+"
          }
          matchMinus = matchPattern(pattern = revCompString(regions.data.frame$consensus.primary[i]), subject = seqC, max.mismatch = maxMis, with.indels = TRUE)
          tempM = as.data.frame(matchMinus)
          tempM$start  = start(matchMinus)
          tempM$end = end(matchMinus)
          if(nrow(tempM) > 0)
          {
            tempM$strand = "-"
          }
          match = rbind(tempP, tempM)
          
          ########## handle overlaps
          primary.size = nchar(regions.data.frame$consensus.primary[i])
          if(nrow(match) > 0)
          {
            match = match[order(match$start),]
            iii = nrow(match)
            
            while(iii > 1)
            {
              overlap = match$start[iii] - match$end[iii - 1] - 1
              
              if(overlap > 0&&overlap <= primary.size%/%10)    #if ovarlap is positive, it means its a gap and small gap will be bridged by adding up to 10 nt to the end of the previous seq
              {
                match$end[iii-1] = match$end[iii-1] + overlap
                
              } else if(overlap < 0 & overlap >= -primary.size) # if overlap is small, seq will be divided evenly between two repeats
              {
                overlap = -overlap
                b = overlap %/% 2
                c = overlap - b
                match$start[iii] = match$start[iii] + b
                match$end[iii - 1] = match$end[iii - 1] - c
                
              } else if(overlap < -primary.size) # if overlap is big, the next seq will be removed completely
              {
                match = match[-c(iii),] 
              }
              iii = iii - 1
            }
          }
          
          if(nrow(match) > 0)
          {
            for(ii in 1 : nrow(match))
            {
              match$seq[ii] = substr(seqC, start = match$start[ii], stop = match$end[ii])
              {
                if(match$strand[ii] == "-")
                {
                  match$seq[ii] = revCompString(match$seq[ii])
                }
              }
            }
            match$width = nchar(match$seq)
            match = match[match$end > match$start,]
            match = match[match$width > 0,]
            match$class = regions.data.frame$class[i]
            match$region.name = paste(assemblyName, fasta.name, sep = "_")
            
          }
          
          #align, need more than 1 sequence
          if(nrow(match) > 1)
          {
            N = regions.data.frame$most.freq.value.N[i]
            if((nrow(match) * N) > LIMIT.REPEATS.TO.ALIGN)
            {
              sample.IDs = sample(x = 1:nrow(match), size = round(LIMIT.REPEATS.TO.ALIGN/N, 0))
            } else
            {
              sample.IDs = 1:nrow(match)
            }
            write.csv(x = match, file = paste(paste(temp.folder, "/", assemblyName, "_out", sep = ""), "/temp2", "/Repeats_", regions.data.frame$name[i], "_", regions.data.frame$start[i], "_", regions.data.frame$end[i], ".csv", sep = ""))
            
            write.fasta(sequences = str_split(match$seq, pattern = "")[sample.IDs], names = paste(assemblyName, ".primary.extract.", regions.data.frame$index[i], ".", regions.data.frame$name[i], ".", match$start, sep = "")[sample.IDs], 
                        file.out = paste(paste(temp.folder, "/", assemblyName, "_out", sep = ""), "/temp2/inputs/primary.extract", ".", regions.data.frame$index[i], ".", regions.data.frame$name[i], ".", "fasta", sep = ""))
            
            
            #mafft the primary extraction
            
            input = paste(paste(temp.folder, "/", assemblyName, "_out", sep = ""), "/temp2/inputs/primary.extract", ".", regions.data.frame$index[i], ".", regions.data.frame$name[i], ".", "fasta", sep = "")
            output = paste(paste(temp.folder, "/", assemblyName, "_out", sep = ""), "/temp2/outputs/primary.extract", ".", regions.data.frame$index[i], ".", regions.data.frame$name[i], ".", "aligned.fasta", sep = "")
            
            if(file.size(input) != 0)
            {
              if(Sys.info()['sysname'] == "Linux")
              {
                system(paste(mafft.bat.file, " --quiet --retree 2 --inputorder ", input, " > ", output, sep = ""), intern = TRUE)
              } 
              
              alignment = read.alignment(output, format = "FASTA", forceToLower = FALSE)
              
              N = regions.data.frame$most.freq.value.N[i]
              if(nchar(alignment$seq[[1]]) < N)
              {
                consensusA = consensus_N(alignment, nchar(alignment$seq[[1]]))
              }
              else
              {
                consensusA = consensus_N(alignment, N)
              }
              
              
              repeats.identified[i] = length(alignment$seq)
              
              consensus.secondary[i] = consensusA
              remove(consensusA)
              
              if(consensus.secondary[i] == "")
              {
                consensus.secondary[i] = "none_identified"
                repeats.identified[i] = 0
              }
            } else
            {
              consensus.secondary[i] = "none_identified"
              repeats.identified[i] = 0
            }
          } 
        } else
        {
          write(paste("on ",regions.data.frame$fasta.name[i], ", skipping consensus for region ", i, "/", nrow(regions.data.frame), " window size: ", regions.data.frame$end[i] - regions.data.frame$start[i], sep = ""), 
                file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
          
          consensus.secondary[i] = "none_identified"
          repeats.identified[i] = 0
        }
      } else
      {
        write(paste("on ",regions.data.frame$fasta.name[i], ", skipping consensus for region ", i, "/", nrow(regions.data.frame), " window size: ", regions.data.frame$end[i] - regions.data.frame$start[i], sep = ""), 
              file = paste(paste(assemblyName, "_out", sep = ""), "/", fasta.name, ".out.txt", sep = ""), append = TRUE)
        
        consensus.secondary[i] = "none_identified"
        repeats.identified[i] = 0
      }
    }
  }
  
  
  if(nrow(regions.data.frame) > 0)
  {
    regions.data.frame = cbind(regions.data.frame, consensus.secondary, repeats.identified)
  }
  write.csv(regions.data.frame, file = paste(temp.folder, "/", assemblyName, "_out/temp1/Temp_late2_regions_", paste(assemblyName, "_out", sep = ""), "_", fasta.name, ".csv", sep = ""))
  
  
  gc()
  
  return(regions.data.frame)
}
