plot_allele_proportions_pop <- function(samples = samples, Reference = Reference, proportion_cutoff = proportion_cutoff, Proportion_variants = Proportion_variants, output_directory = "./allele_proportion_plots"){
  
  # Create a directory for ALL plots
  dir.create(output_directory, showWarnings = FALSE)
  
  Levels_ <- data.frame(rbind(c('Del', 'A', 'C', 'T', 'G'), c(2, 4, 5, 6, 7)),stringsAsFactors=FALSE)
  prop_ALT <- as.data.frame(matrix(NA, nrow=dim(Reference)[1]))
  
  for(i in 1: dim(Reference)[1]){ 
    what_level <- which(Levels_[1,] == Reference$ref[i])
    Levels_temp <- Levels_[,-what_level]
    Prop_temp <- Proportion_variants[i, as.numeric(Levels_temp[2,])]
    prop_ALT[i,1] <- max(Prop_temp)/length(samples)
  }
  
  # Distribute proportions in different columns according to feature coordinates
  # ETS1
  prop_ALT[1:2194,2] <- prop_ALT[1:2194,1]
  # 18S
  prop_ALT[2195:4002,3] <- prop_ALT[2195:4002,1]
  # ITS1
  prop_ALT[4003:4270,4] <- prop_ALT[4003:4270,1]
  # 5.8S
  prop_ALT[4271:4434,5] <- prop_ALT[4271:4434,1]
  # ITS2
  prop_ALT[4435:4622,6] <- prop_ALT[4435:4622,1]
  # 25S
  prop_ALT[4623:8009,7] <- prop_ALT[4623:8009,1]
  # ETS2
  prop_ALT[8010:9378,8] <- prop_ALT[8010:9378,1]
  # Add insertions separately
  prop_ALT[,9] <- Proportion_variants$Ins/length(samples)
  
  pdf(paste(output_directory, '/45SrRNAgene_ALTproportions.pdf', sep=''), width = 24, height = 10)
  par(mar = c(5, 6, 4, 1))
  plot(prop_ALT[,1], cex=0.0, ylim=c(-0.05,1), xlim=c(1,dim(prop_ALT)[1]+700), main = '45S rRNA gene polymorphisms in the population', ylab=paste('Proportion of strains with an alternative allele > ', proportion_cutoff, sep=''), cex.lab=2, cex.axis=2, xlab='Position along the 45S rRNA gene', cex.main=1.8, bty="n")
  lines(prop_ALT[,9], cex=0.1, col='#525252', lty=1, lwd=2)
  lines(prop_ALT[,2], cex=0.1, col='#bdbdbd', lwd=1)
  lines(prop_ALT[,3], cex=0.1, col='#d7191c', lwd=3)
  lines(prop_ALT[,4], cex=0.1, col='#abdda4', lwd=2.5)
  lines(prop_ALT[,5], cex=0.1, col='#fdae61', lwd=3)
  lines(prop_ALT[,6], cex=0.1, col='#abdda4', lwd=2.5)
  lines(prop_ALT[,7], cex=0.1, col='#2b83ba', lwd=3)
  lines(prop_ALT[,8], cex=0.1, col='#bdbdbd')
  rect(2195, -0.08, 4002, 0, col='#d7191c25', border=NA)
  rect(4271, -0.08, 4434, 0, col='#fdae6140', border=NA)
  rect(4623, -0.08, 8009, 0, col='#2b83ba30', border=NA)
  legend("right",  legend=c("18S","5.8S", "25S", "ITS", "ETS", "Insertions"), lty=c(1,1,1,1,1,1), cex=1.5, box.lty=0, text.font=2, text.col =c("#d7191c", "#fdae61", "#2b83ba", "#abdda4", "#bdbdbd", "#525252"), col =c("#d7191c", "#fdae61", "#2b83ba", "#abdda4", "#bdbdbd", "#525252"))
  dev.off()
}