library(dplyr)
library(Biostrings)
library(circlize)
library(tibble)

setwd("~/rds/rds-8b3VcZwY7rY/projects/centromeres/TRASH_final")

species_list <- read.table("genomes.txt")$V1

for(species in species_list){

  name <- species
  
  # Import repeats ----------------------------------------------------------
  
  file <- list.files(paste("runs/", name, sep = ""), pattern = "Summary")
  
  repeats <- read.csv(paste("runs/", name, "/", file, sep = ""))
  
  repeats$most.freq.value.N <- as.integer(repeats$most.freq.value.N)
  
  # Import genome -----------------------------------------------------------
  
  fasta <- list.files(paste("runs/", name, sep = ""), pattern = "\\.fa$")
  
  genome <- readDNAStringSet(paste("runs/", name, "/", fasta, sep = ""))
  
  # Identify Chromosomes ----------------------------------------------------
  
  if(sum(grepl("chr", genome@ranges@NAMES)) > 1){
    chrs <- grep("chr", genome@ranges@NAMES)
  } else if(sum(grepl("SUPER", genome@ranges@NAMES)) > 1){
    temp1 <- grep("SUPER", genome@ranges@NAMES)
    temp2 <- grep("u", genome@ranges@NAMES)
    chrs <- setdiff(temp1, temp2)
  } else {
    chrs <- which(genome@ranges@width > 2500000)
  }
  
  # Identify main peaks -----------------------------------------------------
  
  repeats_l <- filter(repeats, most.freq.value.N > 50)
  
  repeat_no <- aggregate(repeats.identified ~ most.freq.value.N, repeats_l, sum) 
  
  total_length <- repeat_no$most.freq.value.N * repeat_no$repeats.identified
  
  percent_fraction <- total_length*100 / sum(genome@ranges@width)
  
  repeat_sum <- cbind(repeat_no, total_length, percent_fraction) %>% 
    dplyr::rename(., repeat_length = most.freq.value.N) %>% 
    dplyr::rename(., no_repeats = repeats.identified)
  
  peaks <- filter(repeat_sum, percent_fraction >= min(length(chrs) * 0.01, 0.2))
  
  if(nrow(peaks) == 0){
    print(paste("No tandem repeats in", name))
    next
  }
  
  png(file = paste("circos_plots/", name, "/", name, "_peaks.png", sep = ""),
      width = 800, height = 600)
  
  plot(repeat_sum$repeat_length, repeat_sum$no_repeats, type = "h", xlim = c(50, max(peaks$repeat_length)*2),
       xlab = "monomer length", ylab = "no of repeats")
  
  dev.off()
  
  # Merge Peaks -------------------------------------------------------------
  
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
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  
  repeat_merge <- repeat_sum
  
  for(i in 1:nrow(peaks_merge)){
    repeat_merge <- subset(repeat_merge, repeat_length != peaks_merge[i,1] 
                           & repeat_length != peaks_merge[i,1] - 2 & repeat_length != peaks_merge[i,1] - 1
                           & repeat_length != peaks_merge[i,1] + 1 & repeat_length != peaks_merge[i,1] + 2)
  } 
  
  repeat_merge <- rbind(repeat_merge, peaks_merge)
  
  ### export data
  
  name_col <- rep(name, nrow(peaks_merge))
  
  export_peaks <- cbind(name_col, peaks_merge) %>% 
    dplyr::rename(., genome = name_col)
  
  export_peaks <- export_peaks[order(export_peaks$percent_fraction, decreasing = TRUE),]
  
  rownames(export_peaks) <- NULL
  
  write.csv(export_peaks, file = paste("circos_plots/", name, "/", name, "_peaks.csv", sep = ""))
  
  png(file = paste("circos_plots/", name, "/", name, "_peaks_m.png", sep = ""),
      width = 800, height = 600)
  
  plot(repeat_merge$repeat_length, repeat_merge$no_repeats, type = "h", xlim = c(50, max(peaks$repeat_length)*2),
       xlab = "monomer length", ylab = "no of repeats")
  
  dev.off()
  
  # Circos Plots ------------------------------------------------------------
  
  no_peaks <- min(nrow(export_peaks), 10)
  
  circos_repeats <- export_peaks[,2]
  
  ### set sectors ###
  
  repeats_chr <- subset(repeats_l, name %in% genome@ranges@NAMES[chrs])
  names <- genome@ranges@NAMES[chrs]
  start_sect <- rep(0, length(chrs))
  end_sect <- genome@ranges@width[chrs]
  sectors <- data.frame(start_sect, end_sect)
  rownames(sectors) <- names
  
  ##
  
  colours <- c("red", "blue", "green", "purple", "orange", "yellow", "pink", "cyan", "brown", "grey")
  
  ##
    
  ### generate circos plot
  
  cairo_pdf(file = paste("circos_plots/", name, "/", name, "_circos.pdf", sep = ""),
      width = 15, height = 15)
  
  circos.par("track.height" = 0.1, cell.padding = c(0.02, 1.00, 0.02, 1.00))
  circos.initialize(xlim = sectors)
  
  # genome
  circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
    chr=CELL_META$sector.index
    xlim=CELL_META$xlim
    ylim=CELL_META$ylim
    circos.text(mean(xlim), mean(ylim), chr, cex = 1.5)
  }) +
    title(name, cex.main = 2.5, line = -18.5)
  
  circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    circos.axis(h="top")
  })
  
  circos.par("track.height" = 0.15, cell.padding = c(0, 0, 0, 0))
  
  repeat_data <- list()
  
  if(no_peaks < 10){
    for(i in 1:no_peaks){
      repeats <- subset(repeats_chr, most.freq.value.N %in% c(circos_repeats[i], 
                                                              circos_repeats[i] - 2, circos_repeats[i] - 1,
                                                              circos_repeats[i] + 1, circos_repeats[i] + 2)) %>% 
        select(., c(3, 4, 5)) %>% 
        dplyr::rename(., chr = name)
      repeats <- remove_rownames(repeats)
      
      value <- rep(5, nrow(repeats))
      
      repeats <- data.frame(repeats, value)
      
      repeat_data[[i]] <- repeats
    }
  } else {
    for(i in 1:9){
      repeats <- subset(repeats_chr, most.freq.value.N %in% c(circos_repeats[i], 
                                                              circos_repeats[i] - 2, circos_repeats[i] - 1,
                                                              circos_repeats[i] + 1, circos_repeats[i] + 2)) %>% 
        select(., c(3, 4, 5)) %>% 
        dplyr::rename(., chr = name)
      repeats <- remove_rownames(repeats)
      
      value <- rep(5, nrow(repeats))
      
      repeats <- data.frame(repeats, value)
      
      repeat_data[[i]] <- repeats
    }
    i = 10
    repeats <- subset(repeats_chr, most.freq.value.N %in% c(circos_repeats[i:length(circos_repeats)], 
                                                            circos_repeats[i:length(circos_repeats)] - 2, circos_repeats[i:length(circos_repeats)] - 1,
                                                            circos_repeats[i:length(circos_repeats)] + 1, circos_repeats[i:length(circos_repeats)] + 2)) %>% 
      select(., c(3, 4, 5)) %>% 
      dplyr::rename(., chr = name)
    repeats <- remove_rownames(repeats)
    
    value <- rep(5, nrow(repeats))
    
    repeats <- data.frame(repeats, value)
    
    repeat_data[[i]] <- repeats
  }
  
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
  
  # Leftover Repeats --------------------------------------------------------
  
  ### set sectors ###
  
  scaffolds <- setdiff(c(1:length(genome)), chrs)
  
  if(length(scaffolds) == 0){
    dev.off()
    print(paste(name, "done"))
    next
  }
  
  repeats_chr <- subset(repeats_l, name %in% genome@ranges@NAMES[scaffolds])
  names <- genome@ranges@NAMES[scaffolds]
  start_sect <- rep(0, length(scaffolds))
  end_sect <- genome@ranges@width[scaffolds]
  sectors <- data.frame(start_sect, end_sect)
  rownames(sectors) <- names
  
  ##
  
  leftovers <- round(sum(genome@ranges@width[scaffolds])/1000000)
    
  ##
  
  ### generate circos plot
  
  circos.clear()
  
  par(fig = c(0.2, 0.7, 0.25, 0.75), new = T)
  
  tryCatch({
    circos.par("track.height" = 0.02, cell.padding = c(0, 0, 0, 0), gap.degree = 0.5)
    circos.initialize(xlim = sectors)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  # genome
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
                                                            circos_repeats[i] + 1, circos_repeats[i] + 2)) %>% 
      select(., c(3, 4, 5)) %>% 
      dplyr::rename(., chr = name)
    repeats <- remove_rownames(repeats)
    
    value <- rep(5, nrow(repeats))
    
    repeats <- data.frame(repeats, value)
    
    repeat_data[[i]] <- repeats
  }
  
  circos.genomicTrack(repeat_data, ylim=c(0, 1), panel.fun = function(region, value, ...) {
    i = getI(...)
    circos.genomicRect(region, value, ytop = 1, ybottom = 0, col = colours[i], border = NA, ...)
  })
  
  dev.off()
  
  print(paste(name, "done"))
}
