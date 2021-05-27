plot_allele_frequency_indiv <- function(samples = samples, Reference = Reference, output_directory = "./allele_frequency_plots"){
  
  # Create a directory for ALL plots
  dir.create(output_directory, showWarnings = FALSE)
  
  # Loop over each sample
  for(i in 1:length(samples)){
    # Subset complete data set for each sample
    ref <- Reference[ , c(1, grep(pattern = samples[i], x = colnames(Reference))) ]
    
    Levels_ <- data.frame(rbind(c('Del', 'A', 'C', 'T', 'G'), c(3, 5, 6, 7, 8)),stringsAsFactors=FALSE)
    freq_ALT <- as.data.frame(matrix(NA, nrow=dim(ref)[1]))
    
    # Loop to add all alternative alleles for each position
    for(j in 1: dim(ref)[1]){ 
      what_level <- which(Levels_[1,] == ref$ref[j])
      Levels_temp <- Levels_[,-what_level]
      freq_temp <- ref[j, as.numeric(Levels_temp[2,]) ]
      freq_ALT[j,1] <- sum(freq_temp, na.rm = TRUE)
    }
    
    # Distribute frequencies in different columns according to feature coordinates
    # ETS1
    freq_ALT[1:2194,2] <- freq_ALT[1:2194,1]
    # 18S
    freq_ALT[2195:4002,3] <- freq_ALT[2195:4002,1]
    # ITS1
    freq_ALT[4003:4270,4] <- freq_ALT[4003:4270,1]
    # 5.8S
    freq_ALT[4271:4434,5] <- freq_ALT[4271:4434,1]
    # ITS2
    freq_ALT[4435:4622,6] <- freq_ALT[4435:4622,1]
    # 25S
    freq_ALT[4623:8009,7] <- freq_ALT[4623:8009,1]
    # ETS2
    freq_ALT[8010:9378,8] <- freq_ALT[8010:9378,1]
    # Add insertions separately
    freq_ALT[,9] <- ref[ , 4]
    
    pdf(paste(output_directory, '/45SrRNAgene_ALTfreq_', samples[i], '.pdf', sep=''), width = 24, height = 10)
    par(mar = c(5, 6, 4, 1))
    plot(freq_ALT[,1], cex=0.0, ylim=c(-0.05,1), xlim=c(1,dim(freq_ALT)[1]+700), main = paste('45S rRNA gene polymorphisms in ', samples[i], sep=''), ylab='Alternative alleles frequency', cex.lab=2, cex.axis=2, xlab='Position along the 45S rRNA gene', cex.main=1.8, bty="n")
    lines(freq_ALT[,9], cex=0.1, col='#525252', lty=1, lwd=2)
    lines(freq_ALT[,2], cex=0.1, col='#bdbdbd', lwd=1)
    lines(freq_ALT[,3], cex=0.1, col='#d7191c', lwd=3)
    lines(freq_ALT[,4], cex=0.1, col='#abdda4', lwd=2.5)
    lines(freq_ALT[,5], cex=0.1, col='#fdae61', lwd=3)
    lines(freq_ALT[,6], cex=0.1, col='#abdda4', lwd=2.5)
    lines(freq_ALT[,7], cex=0.1, col='#2b83ba', lwd=3)
    lines(freq_ALT[,8], cex=0.1, col='#bdbdbd')
    rect(2195, -0.08, 4002, 0, col='#d7191c25', border=NA)
    rect(4271, -0.08, 4434, 0, col='#fdae6140', border=NA)
    rect(4623, -0.08, 8009, 0, col='#2b83ba30', border=NA)
    abline(h=0.02, lty=2, col="black")
    abline(h=0.05, lty=3, col="black")
    legend("right",  legend=c("18S","5.8S", "25S", "ITS", "ETS", "Insertions", "", "0.05", "0.02"), lty=c(1,1,1,1,1,1,0,3,2), cex=1.5, box.lty=0, text.font=2, text.col =c("#d7191c", "#fdae61", "#2b83ba", "#abdda4", "#bdbdbd", "#525252", NA,  "black", "black"), col =c("#d7191c", "#fdae61", "#2b83ba", "#abdda4", "#bdbdbd", "#525252", NA, "black", "black"))
    dev.off()
  }
}
