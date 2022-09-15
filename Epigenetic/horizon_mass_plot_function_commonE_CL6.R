plot_horizon_mass <- function(cov,patient,peak,smooth,lwd=4,type="histogram",legend){
  combos <- combos[grepl(patient,combos)]
  combos <- combos[grepl("_A|_B|_C|_D|C542_F",combos)]
  regions <- unlist(lapply(1: length(combos), function(c){unlist(str_split(combos[c],"_"))[2]}))
  regions <- c(regions,"E")
  col <- col[regions]
  
  peaks_2000 <- resize(peaks, width = 2000, fix = "center")
  peak_gr <- peaks_2000[peaks_2000$peak == peak,]
  start=start(peaks_2000[peaks_2000$peak == peak,])
  end=end(peaks_2000[peaks_2000$peak == peak,])
  chr=as.character(as.data.frame(peaks_2000)[peaks_2000$peak == peak,1])
  
  #Get the order of presentation so that max value is plotted first
  cov <- c(cov[grepl(patient,names(cov))],cov["pan_E"])
  cov <- cov[grepl("_A|_B|_C|_D|C542_F|pan_E",names(cov))]
  cov<-lapply(cov,subsetByOverlaps,peak_gr)
  suppressWarnings(max<-unlist(lapply(1:length(combos),function(n){
    c <- regions[n]
    max(cov[[c]]$score)})))
  max[which(!is.finite(max))] <- 0
  names(max) <- regions
  max_order <- order(max, decreasing=T)
  max_order <- names(max[max_order])
  
  
  
  # get tracks from coverage and colour them
  tracks <- lapply(1:length(regions), function(c) {
    region <- unlist(str_split(regions[c],"_"))[2]
    mycov <- cov[[c]]
    track <- DataTrack(genome="hg38",
                       range=mycov,
                       fill.mountain="transparent",
                       col.mountain=col,
                       name=regions[c],fill=col,
                       col.histogram=col,
                       fill.histogram="transparent",
                       fill="transparent",
                       col=col,
                       col.line=col,
                       from=start,
                       to = end,
                       chr=chr
    )
    displayPars(track) <- list(groups = factor(region, 
                                               levels = regions),
                               col=col[region],
                               legend = TRUE,
                               missingAsZero=F,
                               lwd=3)
    return(track)})
  
  names(tracks) <- regions
  
  tracks <- tracks[max_order]
  
  names(tracks[[1]]) <- paste0(as.character(patient)," normalised coverage")
  
  ylims <- range(lapply(tracks,values))
  
  tracks_to_plot<-OverlayTrack(tracks,name=as.character(patient),
                               ylim=ylims)
  
  annotracks <- peak_tracks[combos]
  
  logfc <-round(fc[peak,paste0(patient,".pure")],digits = 2)
  sig2 <- sig_filter[peak,paste0(patient,".pure")]
  status<-get_status(patient=patient,
                     peak=peakss$peak[t])
  
  results_plot <- results[results$peak==peak, grepl(patient,colnames(results))]
  results_plot$CvN_pvalue <- p[peak,paste0(patient,".pure")]
  results_plot$CvN_padj <-  p_adj[peak,paste0(patient,".pure")]
  results_plot<-as.data.frame(sapply(results_plot, as.numeric))
  mynames <- rownames(results_plot)
  sig<-ifelse(is.na(results_plot[,1]),"",
              ifelse(results_plot[,1] <= 0.0001 , " ****",
                     ifelse(results_plot[,1] <= 0.001 , " ***",
                            ifelse(results_plot[,1] <= 0.01 , " **",
                                   ifelse(results_plot[,1] <= 0.05 , " *", "")))))
  
  results_plot <- as.data.frame(sapply(sapply(results_plot, as.numeric), scientific))
  results_plot[,1] <- paste0(results_plot[,1],sig)
  rownames(results_plot) <- mynames
  P. <- c(results_plot[c(1,3,5,7),1])
  `Adjusted P.` <- c(results_plot[c(2,4,6,8),1])
  results_plot <- cbind(P.,`Adjusted P.`)
  rownames(results_plot) <- c("~region vs. ~","~purity vs ~","~purity+region vs ~purity","Carcinoma vs. normal")
  
  if(!is.na(sig2) & as.character(sig2)=="TRUE" & get_status(patient,peak)=="Clonal" ){
  plotTracks(background.title="white",fontcolor.title="black",col.axis="black",
             showTitle=T, clip=FALSE,
             c(AnnotationTrack(peaks, name = NULL, fill="grey"),
               annotracks,
               tracks_to_plot),
             legend=T,
             type="l",
             from=start,
             to = end,
             min.distance=10,
             chromosome=chr,
             lwd.baseline=1,
             col.baseline="black",
             windowSize =40,window=-1,
             main=paste0(mcols(peaks)[mcols(peaks)$peak==peak,"symbol"], " ",
                         peak,
                         "\n","Putatively ",get_status(patient,peak),", ",
                        "PASS=",as.character(sig2),", ",
                        "LogFC=",round(fc[peak,paste0(patient,".pure")],digits = 2),
                        ", ",
                        "~purity+region Padj=", results_plot[3,2]
                        ),
             ylim=ylims,cex.title=1,cex.axis=0.7, cex.main=1, add=TRUE,
             sizes=c(0.02,rep(0.02,length(combos)),1),
             detail="coverage")
}

}
