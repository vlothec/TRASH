

circos_plotter <- function(outputs.directory, assemblyName, chr_names, chr_lengths, plot.min, plot.max){
  print("plotting")
  setwd(outputs.directory)
  #file <- list.files(outputs.directory, pattern = "Summary")
  repeats <- read.csv(paste(outputs.directory, "/Summary.of.repetitive.regions.", assemblyName, ".csv", sep = ""))
  repeats$most.freq.value.N <- as.integer(repeats$most.freq.value.N)
  
  if(!dir.exists(paste(outputs.directory, "/plots", sep = "")))
  {
    dir.create(paste(outputs.directory, "/plots", sep = ""))
  }
  print("identify chromosomes")
  ### identify chromosomes
  
  if(sum(grepl("chr", chr_names, ignore.case = TRUE)) > 1){
    
    temp1 <- grep("chr", chr_names, ignore.case = TRUE)
    
    temp2 <- which(chr_lengths > 1000000)
    
    chrs <- intersect(temp1, temp2)
    
  } else if(sum(grepl("SUPER", chr_names, ignore.case = TRUE)) > 1){
    
    temp1 <- grep("SUPER", chr_names, ignore.case = TRUE)
    
    temp2 <- which(chr_lengths > 1000000)
    
    chrs <- intersect(temp1, temp2)
    
  } else {
    
    chrs <- which(chr_lengths > 2500000)
    
  }
  print("Identify main peaks")
  ### Identify main peaks
  
  repeats_l <- repeats[repeats$most.freq.value.N >= plot.min,]
  repeats_l <- repeats[repeats$most.freq.value.N <= plot.max,]
  repeat_no <- aggregate(repeats.identified ~ most.freq.value.N, repeats_l, sum)
  total_length <- repeat_no$most.freq.value.N * repeat_no$repeats.identified
  percent_fraction <- total_length*100 / sum(chr_lengths)
  repeat_sum <- cbind(repeat_no, total_length, percent_fraction)
  colnames(repeat_sum)[1] <- "repeat_length"
  colnames(repeat_sum)[2] <- "no_repeats"
  
  peaks <- repeat_sum[repeat_sum$percent_fraction >= min(length(chrs) * 0.01, 0.2),]
  
  if(nrow(peaks) == 0){
    print(paste("No tandem repeats in", assemblyName))
    return(0)
  }
  
  png(file = paste(outputs.directory, "/plots/", assemblyName, "_peaks.png", sep = ""),
      width = 800, height = 600)
  plot(repeat_sum$repeat_length, repeat_sum$no_repeats, type = "h", xlim = c(50, min(max(peaks$repeat_length)*2, plot.max)),
       xlab = "monomer length", ylab = "no of repeats")
  dev.off()
  
  print("Merge peaks")
  ### Merge peaks
  
  peaks_merge <- peaks[order(peaks$no_repeats, decreasing = TRUE),]
  
  i = 1
  
  tryCatch({
    while(peaks_merge$repeat_length[i] > 0){
      main_peak <- peaks_merge$repeat_length[i]
      j = 1
      for(peak in repeat_sum$repeat_length){
        if(peak == main_peak - 2 | peak == main_peak - 1 | peak == main_peak + 1 | peak == main_peak + 2){
          peaks_merge[i,2] <- peaks_merge[i,2] + repeat_sum[j,2]
          peaks_merge[i,3] <- peaks_merge[i,3] + repeat_sum[j,3]
          peaks_merge[i,4] <- peaks_merge[i,4] + repeat_sum[j,4]
        }
        j = j + 1
      }
      peaks_merge <- subset(peaks_merge, repeat_length != main_peak - 2 & repeat_length != main_peak - 1 &
                              repeat_length != main_peak + 1 & repeat_length != main_peak + 2)
      i = i + 1
    }
  }, error=function(e){})
  
  repeat_merge <- repeat_sum
  
  for(i in 1:nrow(peaks_merge)){
    repeat_merge <- subset(repeat_merge, repeat_length != peaks_merge[i,1] 
                           & repeat_length != peaks_merge[i,1] - 2 & repeat_length != peaks_merge[i,1] - 1
                           & repeat_length != peaks_merge[i,1] + 1 & repeat_length != peaks_merge[i,1] + 2)
  } 
  
  repeat_merge <- rbind(repeat_merge, peaks_merge)
  
  print("export data")
  ### export data
  
  name_col <- rep(assemblyName, nrow(peaks_merge))
  export_peaks <- cbind(name_col, peaks_merge)
  colnames(export_peaks)[1] <- "genome"
  export_peaks <- export_peaks[order(export_peaks$percent_fraction, decreasing = TRUE),]
  rownames(export_peaks) <- NULL
  
  write.csv(export_peaks, file = paste(outputs.directory, "/plots/", assemblyName, "_peaks.csv", sep = ""))
  png(file = paste(outputs.directory, "/plots/", assemblyName, "_peaks_m.png", sep = ""),
      width = 800, height = 600)
  plot(repeat_merge$repeat_length, repeat_merge$no_repeats, type = "h", xlim = c(50, min(max(peaks$repeat_length)*2, 2000)),
       xlab = "monomer length", ylab = "no of repeats")
  dev.off()
  
  print("Circos Plotting")
  ### Circos Plotting
  
  no_peaks <- min(nrow(export_peaks), 10)
  circos_repeats <- export_peaks[,2]
  colours <- c("red", "blue", "green", "purple", "orange", "yellow", "pink", "cyan", "brown", "grey")
  
  print("set sectors")
  # set sectors
  
  repeats_chr <- subset(repeats_l, name %in% chr_names[chrs])
  names <- chr_names[chrs]
  start_sect <- rep(0, length(chrs))
  end_sect <- chr_lengths[chrs]
  sectors <- data.frame(start_sect, end_sect)
  rownames(sectors) <- names
  
  print("generate circos plot")
  # generate circos plot
  
  cairo_pdf(file = paste(outputs.directory, "/plots/", assemblyName, "_circos.pdf", sep = ""),
            width = 15, height = 15)
  
  print("1")
  circos.par("track.height" = 0.1, cell.padding = c(0.02, 1.00, 0.02, 1.00))
  circos.initialize(xlim = sectors)
  
  print("2")
  circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
    chr=CELL_META$sector.index
    xlim=CELL_META$xlim
    ylim=CELL_META$ylim
    circos.text(mean(xlim), mean(ylim), chr, cex = 1.5)
  }) +
    title(assemblyName, cex.main = 2.5, line = -18.5)
  
  print("3")
  circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    circos.axis(h="top")
  })
  
  print("4")
  circos.par("track.height" = 0.15, cell.padding = c(0, 0, 0, 0))
  
  repeat_data <- list()
  
  print("5")
  if(no_peaks < 10){
    for(i in 1:no_peaks){
      repeats <- subset(repeats_chr, most.freq.value.N %in% c(circos_repeats[i], 
                                                              circos_repeats[i] - 2, circos_repeats[i] - 1,
                                                              circos_repeats[i] + 1, circos_repeats[i] + 2))[,c(2,3,4)]  
      colnames(repeats)[1] <- "chr"
      row.names(repeats) <- NULL
      value <- rep(5, nrow(repeats))
      repeats <- data.frame(repeats, value)
      repeat_data[[i]] <- repeats
    }
  } else {
    for(i in 1:9){
      repeats <- subset(repeats_chr, most.freq.value.N %in% c(circos_repeats[i], 
                                                              circos_repeats[i] - 2, circos_repeats[i] - 1,
                                                              circos_repeats[i] + 1, circos_repeats[i] + 2))[,c(2,3,4)]  
      colnames(repeats)[1] <- "chr"
      row.names(repeats) <- NULL
      value <- rep(5, nrow(repeats))
      repeats <- data.frame(repeats, value)
      repeat_data[[i]] <- repeats
    }
    i = 10
    repeats <- subset(repeats_chr, most.freq.value.N %in% c(circos_repeats[i], 
                                                            circos_repeats[i] - 2, circos_repeats[i] - 1,
                                                            circos_repeats[i] + 1, circos_repeats[i] + 2))[,c(2,3,4)]  
    colnames(repeats)[1] <- "chr"
    row.names(repeats) <- NULL
    value <- rep(5, nrow(repeats))
    repeats <- data.frame(repeats, value)
    repeat_data[[i]] <- repeats
  }
  
  print("genomic track")
  circos.genomicTrack(repeat_data, ylim=c(0, 1), panel.fun = function(region, value, ...) {
    i = getI(...)
    circos.genomicRect(region, value, ytop = 1, ybottom = 0, col = colours[i], border = NA, ...)
  })
  
  if(no_peaks < 10){
    legend(x = 0.425, y = 0, legend = c(export_peaks[1:no_peaks,2]),
           col = colours[1:no_peaks], lty=1, lwd = 2.5, yjust = 0.5, cex = 1.5)
  } else {
    legend(x = 0.425, y = 0, legend = c(export_peaks[1:(no_peaks-1),2], "rest"),
           col = colours[1:10], lty=1, lwd = 2.5, yjust = 0.5, cex = 1.5)
  }
  
  ### leftover repeats
  print("leftover repeats")
  
  scaffolds <- setdiff(c(1:length(chr_names)), chrs)
  if(length(scaffolds) == 0){
    dev.off()
    print(paste(assemblyName, " plots done", sep = ""))
    return(0)
  }
  repeats_chr <- subset(repeats_l, name %in% chr_names[scaffolds])
  names <- chr_names[scaffolds]
  start_sect <- rep(0, length(scaffolds))
  end_sect <- chr_lengths[scaffolds]
  sectors <- data.frame(start_sect, end_sect)
  rownames(sectors) <- names
  
  leftovers <- round(sum(chr_lengths[scaffolds])/1000000)
  
  # plot leftovers
  print("plot leftovers")
  
  circos.clear()
  
  par(fig = c(0.2, 0.7, 0.25, 0.75), new = T)
  
  tryCatch({
    circos.par("track.height" = 0.02, cell.padding = c(0, 0, 0, 0), gap.degree = 0.5)
    circos.initialize(xlim = sectors)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
    chr=CELL_META$sector.index
    xlim=CELL_META$xlim
    ylim=CELL_META$ylim
  }) +
    title(paste("leftover scaffolds = ", leftovers, "Mb", sep = ""), line = -18, cex.main = 1.5)
  
  circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    circos.axis(h="top", labels.cex = 0.5)
  })
  
  circos.par("track.height" = 0.2, cell.padding = c(0, 0, 0, 0))
  
  repeat_data <- list()
  
  for(i in 1:no_peaks){
    repeats <- subset(repeats_chr, most.freq.value.N %in% c(circos_repeats[i], 
                                                            circos_repeats[i] - 2, circos_repeats[i] - 1,
                                                            circos_repeats[i] + 1, circos_repeats[i] + 2))[,c(2,3,4)]  
    colnames(repeats)[1] <- "chr"
    row.names(repeats) <- NULL
    value <- rep(5, nrow(repeats))
    repeats <- data.frame(repeats, value)
    repeat_data[[i]] <- repeats
  }
  
  circos.genomicTrack(repeat_data, ylim=c(0, 1), panel.fun = function(region, value, ...) {
    i = getI(...)
    circos.genomicRect(region, value, ytop = 1, ybottom = 0, col = colours[i], border = NA, ...)
  })
  
  dev.off()
  
  print(paste(assemblyName, " plots done", sep = ""))
  
  return(0)
}