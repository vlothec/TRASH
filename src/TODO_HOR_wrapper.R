HOR.wrapper = function(temp.folder = "", 
                       assemblyName = "", 
                       mafft.bat.file = "",
                       TRESHOLD = "",
                       CUTOFF = "")
{
  setwd(temp.folder)
  
  repeats.files = list.files(path = cen180.files.directory, pattern = "all.repeats.from", recursive = FALSE) 
  
  if(HOR.building.strategy == "individual")
  {
    pairs.to.compare = data.frame(n = vector(mode = "numeric", length = 0), m = vector(mode = "numeric", length = 0))
    l = 1
    for(n in 1 : 5)
    {
      for(m in n : 5)
      {
        pairs.to.compare = rbind(pairs.to.compare, c(n,m))
        l = l + 1
      }
    }
    names(pairs.to.compare) = c("m", "n")
    
    #foreach(i = 1 : length(repeats.files)) %dopar% {
    for(i in 1 : 66){#1
      
      repeats = read.csv(paste(cen180.files.directory, repeats.files[i], sep = "/"), stringsAsFactors = FALSE)
      repeats$class[is.na(repeats$class)] = ""
      repeats = repeats[repeats$class == "cen180",]
      repeats$strand[repeats$strand == "+"] = "1"
      repeats$strand[repeats$strand == "-"] = "2"
      
      
      #foreach(j = 1 : nrow(pairs.to.compare)) %dopar% {
      #for(j in 1 : 15){
      j = as.numeric(task_id)
      print(paste("working on ", j, " pair from ", repeats.files[i], sep = ""))
      
      if(pairs.to.compare$m[j] == pairs.to.compare$n[j])
      {
        repeatsT = repeats[repeats$chromosome == unique(repeats$chromosome)[pairs.to.compare$m[j]],]
        
        file.to.align = paste("to_align_", unique(repeats$chromosome)[pairs.to.compare$m[j]], "_", repeats.files[i], ".fasta", sep = "")
        
        write.fasta(sequences = as.list(repeatsT$seq), names = as.list(paste(1:length(repeatsT$strand), "_D", repeatsT$strand,  sep = "")), file.out = file.to.align, open = "w")
        
        output.alignment = paste("aligned_", unique(repeats$chromosome)[pairs.to.compare$m[j]], "_", repeats.files[i], ".fasta", sep = "")
        
        system(paste("mafft --quiet --inputorder --kimura 1 --retree 1 ", file.to.align, " > ", output.alignment, sep = ""), intern = TRUE)
        
        system(paste("/home/pw457/HOR/HOR.V3.3", " ", temp.folder, "/ ", output.alignment, " ", 1, " ", threshold, " ", cutoff, " ", 1, sep = ""), intern = FALSE)
        
        hors = read.csv(file = paste(temp.folder, "/", "HORs_method_1_", output.alignment, "_t_", threshold, "_c_", cutoff, ".csv", sep = ""), header = TRUE)
        
        if(nrow(hors) > 1)
        {
          hors = as.data.frame(hors)
          
          hors$start.A.bp = repeatsT$start[hors$start_A]
          hors$start.B.bp = repeatsT$start[hors$start_B]
          hors$end.A.bp = repeatsT$end[hors$end_A]
          hors$end.B.bp = repeatsT$end[hors$end_B] 
          hors$chrA = repeatsT$chromosome[hors$start_A]  
          hors$chrB = repeatsT$chromosome[hors$start_B]
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
          
          png(filename = paste(temp.folder, "/HOR_plot_", unique(repeats$chromosome)[pairs.to.compare$n[j]], "_", repeats.files[i], ".png", sep = ""), width = 5000, height = 5000, pointsize = 25)
          plot(xlab = unique(repeats$chromosome)[pairs.to.compare$m[j]], ylab = unique(repeats$chromosome)[pairs.to.compare$n[j]], 
               x = hors$start.A.bp, pch = 19, cex = 0.1, 
               y = hors$start.B.bp, 
               xlim = c(min(hors$start.A.bp), max(hors$start.A.bp)), 
               ylim = c(min(hors$start.B.bp), max(hors$start.B.bp)))
          dev.off()
        }
        
        
        remove(hors, repeatsT)
        gc()
      } else
      {
        time1 = Sys.time()
        repeatsN = repeats[repeats$chromosome == unique(repeats$chromosome)[pairs.to.compare$n[j]],]
        repeatsM = repeats[repeats$chromosome == unique(repeats$chromosome)[pairs.to.compare$m[j]],]
        repeatsT = rbind(repeatsN, repeatsM)
        
        split.no = nrow(repeatsN)
        
        file.to.align = paste("to_align_", unique(repeats$chromosome)[pairs.to.compare$m[j]], "_", unique(repeats$chromosome)[pairs.to.compare$n[j]], "_", repeats.files[i], ".fasta", sep = "")
        
        write.fasta(sequences = as.list(repeatsT$seq), names = as.list(paste(1:length(repeatsT$strand), "_D", repeatsT$strand,  sep = "")), file.out = file.to.align, open = "w")
        
        output.alignment = paste("aligned_", unique(repeats$chromosome)[pairs.to.compare$m[j]], "_", unique(repeats$chromosome)[pairs.to.compare$n[j]], "_", repeats.files[i], ".fasta", sep = "")
        
        system(paste("mafft --quiet --inputorder --kimura 1 --retree 1 ", file.to.align, " > ", output.alignment, sep = ""), intern = TRUE)
        
        system(paste("/home/pw457/HOR/HOR.V3.3", " ", temp.folder, "/ ", output.alignment, " ", split.no, " ", threshold, " ", cutoff, " ", 1, sep = ""), intern = TRUE)
        
        
        hors = read.csv(file = paste(temp.folder, "/", "HORs_method_1_", output.alignment, "_t_", threshold, "_c_", cutoff, ".csv", sep = ""), header = TRUE)
        
        hors = as.data.frame(hors)
        if(nrow(hors) > 1)
        {
          hors$start.A.bp = repeatsT$start[hors$start_A]
          hors$start.B.bp = repeatsT$start[hors$start_B]
          hors$end.A.bp = repeatsT$end[hors$end_A]
          hors$end.B.bp = repeatsT$end[hors$end_B] 
          hors$chrA = repeatsT$chromosome[hors$start_A]  
          hors$chrB = repeatsT$chromosome[hors$start_B]
          hors$block.size.in.units = hors$end_A - hors$start_A + 1
          hors$block.A.size.bp = hors$end.A.bp - hors$start.A.bp + 1
          hors$block.B.size.bp = hors$end.B.bp - hors$start.B.bp + 1
          for(l in nrow(hors) : 2)
          {
            block1size = abs(hors$end.A.bp[l] - hors$start.A.bp[l])
            block2size = abs(hors$end.B.bp[l] - hors$start.B.bp[l])
            if(abs(block2size - block1size) > 1000)
            {
              hors = hors[-l,]
              l = l - 1
            }
          }
          
          
          png(filename = paste(temp.folder, "/HOR_plot_", 
                               unique(repeats$chromosome)[pairs.to.compare$m[j]], "_", 
                               unique(repeats$chromosome)[pairs.to.compare$n[j]], "_", repeats.files[i],
                               length(hors$start.A.bp[hors$start_A < split.no & hors$start_B < split.no])/nrow(hors),
                               ".png", sep = ""), width = 15000, height = 15000, pointsize = 25)
          par(mfrow = c(2,2))
          
          plot(xlab = unique(repeats$chromosome)[pairs.to.compare$m[j]], ylab = unique(repeats$chromosome)[pairs.to.compare$n[j]], 
               x = hors$start.A.bp[hors$start_A < split.no & hors$start_B > split.no], pch = 19, cex = 0.1, 
               y = hors$start.B.bp[hors$start_A < split.no & hors$start_B > split.no], 
               xlim = c(min(hors$start.A.bp), max(hors$start.A.bp)), 
               ylim = c(min(hors$start.B.bp), max(hors$start.B.bp)))
          
          plot(xlab = unique(repeats$chromosome)[pairs.to.compare$m[j]], ylab = unique(repeats$chromosome)[pairs.to.compare$n[j]], 
               x = hors$start.A.bp[hors$start_A < split.no & hors$start_B < split.no], pch = 19, cex = 0.1, 
               y = hors$start.B.bp[hors$start_A < split.no & hors$start_B < split.no], 
               xlim = c(min(hors$start.A.bp), max(hors$start.A.bp)), 
               ylim = c(min(hors$start.B.bp), max(hors$start.B.bp)))
          
          plot(xlab = unique(repeats$chromosome)[pairs.to.compare$m[j]], ylab = unique(repeats$chromosome)[pairs.to.compare$n[j]], 
               x = hors$start.A.bp[hors$start_A > split.no & hors$start_B > split.no], pch = 19, cex = 0.1, 
               y = hors$start.B.bp[hors$start_A > split.no & hors$start_B > split.no], 
               xlim = c(min(hors$start.A.bp), max(hors$start.A.bp)), 
               ylim = c(min(hors$start.B.bp), max(hors$start.B.bp)))
          
          dev.off()
        }
        
        write.csv(x = hors, file = paste("Formatted.HORs_", unique(repeats$chromosome)[pairs.to.compare$m[j]], "_", unique(repeats$chromosome)[pairs.to.compare$n[j]], "_", repeats.files[i], sep = ""))
        remove(hors, repeatsT)
        gc()
        Sys.time() = time1
      }
      
      do.overlaps = FALSE
      if(do.overlaps)
      {
        # check overlapping HORs and remove them
        elo = 0
        for( i in 1 : nrow(hors))
        {
          if(hors$direction.1.para_2.perp.[i] == 2)
          {
            if(hors$start_B[i] - hors$end_A[i] < 1)
            {
              elo = c(elo, i)
              #hors = hors[-i,]
            }
          }
        }
        elo = elo[-1]
        #1568 cases
        if(length(elo) > 0)
        {
          hors = hors[-elo,]
        }
        
        elo = 0 #check HORs that overlap in forward direction
        for(i in 1 : nrow(hors))
        {
          if(hors$direction.1.para_2.perp.[i] == 1)
          {
            if(hors$start_B[i] - hors$end_A[i] < 1)
            {
              elo = c(elo, i)
            }
          }
        }
        elo = elo[-1]
        #0 cases
      }
      
      dump.lines = FALSE
      if(dump.lines)
      {
        
        png(filename = paste(temp.folder, "/HOR_plot_", unique(repeats$chromosome)[pairs.to.compare$m[j]], "_", unique(repeats$chromosome)[pairs.to.compare$n[j]], "_", repeats.files[i], ".png", sep = ""), 
            width = 5000, height = 5000, pointsize = 25)
        plot(xlab = unique(repeats$chromosome)[pairs.to.compare$m[j]], ylab = unique(repeats$chromosome)[pairs.to.compare$n[j]], 
             x = hors$start.A.bp, pch = 19, cex = 0.1, 
             y = hors$start.B.bp, 
             xlim = c(min(hors$start.A.bp), max(hors$start.A.bp)), 
             ylim = c(min(hors$start.B.bp), max(hors$start.B.bp)))
        dev.off()
        
        #count HORs
        hors = read.csv(file = "C:/Users/wlodz/Desktop/cenRepeats/HOR.wrapper/HORs_method_1_aligned_Chr1_Chr3_ed25.11.21.all.repeats.from.t2t-col.20210610.fasta.csv.fasta_t_5_c_3.csv", header = TRUE)
        png(filename = paste("C:/Users/wlodz/Desktop/cenRepeats/HOR.wrapper/HOR_plot_", unique(repeats$chromosome)[pairs.to.compare$m[j]], "_", unique(repeats$chromosome)[pairs.to.compare$n[j]], "_", repeats.files[i],length(hors$start.A.bp[hors$start_A < split.no & hors$start_B < split.no])/nrow(hors), ".png", sep = ""), width = 15000, height = 15000, pointsize = 25)
        
        # make activity (count x size) plots
        
        cen180$HORcount = 0
        cen180$HORlengthsSum = 0
        
        for(i in 1 : length(hors$index))
        {
          print(i)
          for(j in 0 : (hors$block.size.in.units[i] - 1))
          {
            cen180$HORcount[hors$START.A[i] + j] = cen180$HORcount[hors$START.A[i] + j] + 1
            cen180$HORcount[hors$START.B[i] + j] = cen180$HORcount[hors$START.B[i] + j] + 1
            cen180$HORlengthsSum[hors$START.A[i] + j] = cen180$HORlengthsSum[hors$START.A[i] + j] + hors$block.size.in.units[i]
            cen180$HORlengthsSum[hors$START.B[i] + j] = cen180$HORlengthsSum[hors$START.B[i] + j] + hors$block.size.in.units[i]
          }
        }
        
        
        cen180$averagedSUMofHORlengths = 0
        for(i in 25 : (nrow(cen180) - 25))
        {
          cen180$averagedSUMofHORlengths[i] = ave(c(cen180$HORlengthsSum[i - 24]:cen180$HORlengthsSum[i + 24]))[1]
        }
        
        png(filename = "C:/Users/wlodz/Desktop/T2T0610 activity new.png", width = 2000, height = 2500, pointsize = 25)
        par(mfrow=c(5,1))
        plot(cen180$start[cen180$chromosome == 1], cen180$averagedSUMofHORlengths[cen180$chromosome == 1], cex.lab = 1.5, cex.axis = 1.5, cex = 0.5, pch = 19, ylim = c(0,max(cen180$averagedSUMofHORlengths)+100), main ="Chromosome 1 repeats. Sum of HORs sizes each repeat is a part of", ylab = "Activity score", xlab = "Position, bp")
        plot(cen180$start[cen180$chromosome == 2], cen180$averagedSUMofHORlengths[cen180$chromosome == 2], cex.lab = 1.5, cex.axis = 1.5, cex = 0.5, pch = 19, ylim = c(0,max(cen180$averagedSUMofHORlengths)+100),  main ="Chromosome 2 repeats. Sum of HORs sizes each repeat is a part of", ylab = "Activity score", xlab = "Position, bp")
        plot(cen180$start[cen180$chromosome == 3], cen180$averagedSUMofHORlengths[cen180$chromosome == 3], cex.lab = 1.5, cex.axis = 1.5, cex = 0.5, pch = 19, ylim = c(0,max(cen180$averagedSUMofHORlengths)+100),  main ="Chromosome 3 repeats. Sum of HORs sizes each repeat is a part of", ylab = "Activity score", xlab = "Position, bp")
        plot(cen180$start[cen180$chromosome == 4], cen180$averagedSUMofHORlengths[cen180$chromosome == 4], cex.lab = 1.5, cex.axis = 1.5, cex = 0.5, pch = 19, ylim = c(0,max(cen180$averagedSUMofHORlengths)+100),  main ="Chromosome 4 repeats. Sum of HORs sizes each repeat is a part of", ylab = "Activity score", xlab = "Position, bp")
        plot(cen180$start[cen180$chromosome == 5], cen180$averagedSUMofHORlengths[cen180$chromosome == 5], cex.lab = 1.5, cex.axis = 1.5, cex = 0.5, pch = 19, ylim = c(0,max(cen180$averagedSUMofHORlengths)+100),  main ="Chromosome 5 repeats. Sum of HORs sizes each repeat is a part of", ylab = "Activity score", xlab = "Position, bp")
        dev.off()
        
      }
      
      old.HOR.2.to.be.fixed = FALSE
      if(old.HOR.2.to.be.fixed)
      {
        time1 = Sys.time()
        repeatsN = repeats[repeats$chromosome == unique(repeats$chromosome)[pairs.to.compare$n[j]],]
        repeatsM = repeats[repeats$chromosome == unique(repeats$chromosome)[pairs.to.compare$m[j]],]
        repeatsT = rbind(repeatsN, repeatsM)
        
        split.no = nrow(repeatsN)
        
        file.to.align = paste("to_align_", unique(repeats$chromosome)[pairs.to.compare$m[j]], "_", unique(repeats$chromosome)[pairs.to.compare$n[j]], "_", repeats.files[i], ".fasta", sep = "")
        
        write.fasta(sequences = as.list(repeatsT$seq), names = as.list(paste(1:length(repeatsT$strand), "_D", repeatsT$strand,  sep = "")), file.out = file.to.align, open = "w")
        
        output.alignment = paste("aligned_", unique(repeats$chromosome)[pairs.to.compare$m[j]], "_", unique(repeats$chromosome)[pairs.to.compare$n[j]], "_", repeats.files[i], ".fasta", sep = "")
        
        system(paste("mafft --quiet --inputorder --kimura 1 --retree 1 ", file.to.align, " > ", output.alignment, sep = ""), intern = TRUE)
        
        system(paste("/home/pw457/HOR/HOR.V3.3", " ", temp.folder, "/ ", output.alignment, " ", split.no, " ", threshold, " ", cutoff, " ", 2, sep = ""), intern = TRUE)
        
        
        hors = read.csv(file = paste(temp.folder, "/", "HORs_method_2_", output.alignment, "_t_", threshold, "_c_", cutoff, ".csv", sep = ""), header = TRUE)
        
        hors = as.data.frame(hors)
        if(nrow(hors) > 1)
        {
          hors$start.A.bp = repeatsT$start[hors$start_A]
          hors$start.B.bp = repeatsT$start[hors$start_B]
          hors$end.A.bp = repeatsT$end[hors$end_A]
          hors$end.B.bp = repeatsT$end[hors$end_B] 
          hors$chrA = repeatsT$chromosome[hors$start_A]  
          hors$chrB = repeatsT$chromosome[hors$start_B]
          hors$block.size.in.units = hors$end_A - hors$start_A + 1
          hors$block.A.size.bp = hors$end.A.bp - hors$start.A.bp + 1
          hors$block.B.size.bp = hors$end.B.bp - hors$start.B.bp + 1
          for(i in nrow(hors) : 1)
          {
            block1size = abs(hors$end.A.bp[i] - hors$start.A.bp[i])
            block2size = abs(hors$end.B.bp[i] - hors$start.B.bp[i])
            if(abs(block2size - block1size) > 1000)
            {
              print(i)
              hors = hors[-i,]
              i = i - 1
            }
          }
          
          png(filename = paste(temp.folder, "/HOR_plot_", unique(repeats$chromosome)[pairs.to.compare$m[j]], "_", unique(repeats$chromosome)[pairs.to.compare$n[j]], "_", repeats.files[i], ".png", sep = ""), width = 5000, height = 5000, pointsize = 25)
          plot(xlab = unique(repeats$chromosome)[pairs.to.compare$m[j]], ylab = unique(repeats$chromosome)[pairs.to.compare$n[j]], 
               x = hors$start.A.bp, pch = 19, cex = 0.1, 
               y = hors$start.B.bp, 
               xlim = c(min(hors$start.A.bp), max(hors$start.A.bp)), 
               ylim = c(min(hors$start.B.bp), max(hors$start.B.bp)))
          dev.off()
        }
      }
    }  
  }
}