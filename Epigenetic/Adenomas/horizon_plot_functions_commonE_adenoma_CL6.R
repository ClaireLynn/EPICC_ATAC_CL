plot_horizon <- function(cov,tracks,patient,peak,smooth,lwd=4,type="histogram",legend){
  combos <- combos[grepl(patient,combos)]
  regions <- unlist(lapply(1: length(combos), function(c){unlist(str_split(combos[c],"_"))[2]}))
  
  peaks_2000 <- resize(peaks, width = 2000, fix = "center")
  peak_gr <- peaks_2000[peaks_2000$peak == peak,]
  start=start(peaks_2000[peaks_2000$peak == peak,])
  end=end(peaks_2000[peaks_2000$peak == peak,])
  chr=as.character(as.data.frame(peaks_2000)[peaks_2000$peak == peak,1])
  
  #Get the order of presentation so that max value is plotted first
  cov <- cov[grepl(patient,names(cov))]
  cov<-lapply(cov,subsetByOverlaps,peak_gr)
  suppressWarnings(max<-unlist(lapply(1:length(combos),function(n){
    c <- combos[n]
    max(cov[[c]]$score)})))
  max[which(!is.finite(max))] <- 0
  names(max) <- combos
  max_order <- order(max, decreasing=T)
  max_order <- names(max[max_order])
  
  # get tracks from coverage and colour them
  tracks <- lapply(1:length(combos), function(c) {
    region <- unlist(str_split(combos[c],"_"))[2]
    mycov <- cov[[c]]
    track <- DataTrack(range=mycov,
                       fill.mountain="transparent",
                       col.mountain=col,
      name=combos[c],fill=col,
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
  
  names(tracks) <- combos
  
  tracks <- tracks[max_order]
  
  names(tracks[[1]]) <- paste0(as.character(patient)," normalised coverage")
  
  ylims <- range(lapply(tracks,values))
  
  tracks_to_plot<-OverlayTrack(tracks,name=as.character(patient),
                               ylim=ylims)
  
  annotracks <- peak_tracks[combos]
  

  plotTracks(background.title="white",fontcolor.title="black",col.axis="black",
             showTitle=T, clip=FALSE,
             c(AnnotationTrack(peaks, name = NULL, fill="grey"),
               annotracks,
               tracks_to_plot),
             legend=legend,
             type="l",
             from=start,
             to = end,
             min.distance=10,
             chromosome=chr,
             lwd.baseline=1,
             col.baseline="black",
             windowSize =40,window=-1,
             main=paste(mcols(peaks)[mcols(peaks)$peak==peak,"id"],
                        mcols(peaks)[mcols(peaks)$peak==peak,"type"] ,
                        mcols(peaks)[mcols(peaks)$peak==peak,"event_type_chr"],"\n",
                        peak, 
                        sep=" "),
             ylim=ylims,cex.title=1,cex.axis=0.7, cex.main=1, add=TRUE,
             sizes=c(0.02,rep(0.02,length(combos)),1),
             detail="coverage")
}



col <- brewer.pal(name="Set1", n=9)
names(col) <- c("A","B","C","D","E","F","G","H","I")


plot_purity <- function(patient=patient){
  purity_plot <- coldata[coldata$patient==patient & !grepl("E",coldata$region),]
  purity_plot$region <- as.factor(as.character(purity_plot$region)) 
  ggplot(purity_plot, aes(region,purity, fill=region), ylim=c(0,1)) + 
    geom_boxplot() +
    theme_classic() +
    scale_y_continuous(limits = c(0,1.0), expand=c(0,0),labels = scales::percent) +
    theme(axis.line.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "none") +
    scale_fill_manual(values = col)
}



Pval_table <- function(patient=patient,peak=peak){
  colname <- paste0(patient,".adenoma")
  NvA_logfc <- round(fc[peak,colname], digits=2)
  NvA_pvalue <- p[peak,colname]
  NvA_padj <- p_adj[peak,colname]
  subtissue_logfc <- subtissue_nf[subtissue_nf$peak==peak,paste0(patient,"_logfc")]
  subtissue_pvalue <- subtissue_nf[subtissue_nf$peak==peak,paste0(patient,"_pvalue")]
  subtissue_padj <- subtissue_nf[subtissue_nf$peak==peak,paste0(patient,"_padj")]
  results_plot <- rbind(NvA_pvalue,NvA_padj,subtissue_pvalue,subtissue_padj)
  results_plot<-as.data.frame(as.numeric(results_plot[,1]))
  sig<-ifelse(is.na(results_plot[,1]),"",
                ifelse(results_plot[,1] <= 0.0001 , " ****",
                  ifelse(results_plot[,1] <= 0.001 , " ***",
                     ifelse(results_plot[,1] <= 0.01 , " **",
                        ifelse(results_plot[,1] <= 0.05 , " *", "")))))
  
  results_plot <- as.data.frame(sapply(sapply(results_plot, as.numeric), scientific))
  results_plot[,1] <- paste0(results_plot[,1],sig)
  results_plot<-rbind(results_plot,NvA_logfc,subtissue_logfc)
  LogFC <- c(results_plot[c(5,6),1])
  P. <- c(results_plot[c(1,3),1])
  `Adjusted P.` <- c(results_plot[c(2,4),1])
  results_plot <- cbind(LogFC,P.,`Adjusted P.`)
  rownames(results_plot) <- c("Adenoma vs. normal","Adenoma vs. Cancer\n~purity+type")
  
  title <- textGrob(get_status_adenoma(patient,peak),
                    gp=gpar(fontsize=12))
  
  padding <- unit(5,"mm")
  t1<-tableGrob(results_plot, theme=ttheme_default(base_size=10))
  table <- gtable_add_rows(
    t1, 
    heights = grobHeight(title) + padding,
    pos = 0)
  table <- gtable_add_grob(
    table, 
    title, 
    1, 1,1, ncol(table))
  grid.draw(table)
}


get_status_adenoma <- function(patient=patient,peak=peak){
  colname <- paste0(patient,".adenoma")
  subtissue_padj <- subtissue_nf[subtissue_nf$peak==peak,paste0(patient,"_padj")]
  NvA_padj <- p_adj[peak,colname]
  if(is.na(subtissue_padj) & is.na(NvA_padj)){title <- "NA"
  } else if(!is.na(subtissue_padj) & subtissue_padj < 0.05){ title <- "Putatively Subclonal"
  } else if(sig_filter[peak,colname]==TRUE){ title <- "Putatively Clonal"
          } else { title <- "NS"}
  return(title)
}


vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
plot_info <- function(pdf_name, patient, peak, type){
  pdf(pdf_name, height=5, width=10.4)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(3, 2, widths = c(6,4.4))))
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  pushViewport(viewport(layout.pos.row=1:3,layout.pos.col=1))
  print(plot_horizon(cov=cov,
                     patient=patient,
                     peak=peak,
                     smooth=30,
                     lwd=1.5,
                     legend=T,
                     type="l"), vp= vplayout(x=1:3,y=1))
  popViewport()
  pushViewport(viewport(layout.pos.row=1,layout.pos.col=2))
  print(Pval_table(patient=patient,peak=peak),vp = vplayout(3, 2))
  popViewport()
  pushViewport(viewport(layout.pos.row=2:3,layout.pos.col=2))
  p1<-plot_purity(patient=patient)
  print(p1,vp = vplayout(2:3,2))
  popViewport()
  graphics.off()
}



#C561 is a special case with 2 adenomas, so here are specific functions for it: 

plot_horizon_C561 <- function(cov,tracks,patient,peak,smooth,lwd=4,type="histogram",legend){
  combos <- combos[grepl(patient,combos)]
  regions <- unlist(lapply(1: length(combos), function(c){unlist(str_split(combos[c],"_"))[2]}))
  
  peaks_2000 <- resize(peaks, width = 2000, fix = "center")
  peak_gr <- peaks_2000[peaks_2000$peak == peak,]
  start=start(peaks_2000[peaks_2000$peak == peak,])
  end=end(peaks_2000[peaks_2000$peak == peak,])
  chr=as.character(as.data.frame(peaks_2000)[peaks_2000$peak == peak,1])
  
  #Get the order of presentation so that max value is plotted first
  cov <- cov[grepl(patient,names(cov))]
  cov<-lapply(cov,subsetByOverlaps,peak_gr)
  suppressWarnings(max<-unlist(lapply(1:length(combos),function(n){
    c <- combos[n]
    max(cov[[c]]$score)})))
  max[which(!is.finite(max))] <- 0
  names(max) <- combos
  max_order <- order(max, decreasing=T)
  max_order <- names(max[max_order])
  
  # get tracks from coverage and colour them
  tracks <- lapply(1:length(combos), function(c) {
    region <- unlist(str_split(combos[c],"_"))[2]
    mycov <- cov[[c]]
    track <- DataTrack(range=mycov,
                       fill.mountain="transparent",
                       col.mountain=col,
                       name=combos[c],fill=col,
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
  
  names(tracks) <- combos
  
  tracks <- tracks[max_order]
  
  names(tracks[[1]]) <- paste0(as.character(patient)," normalised coverage")
  
  ylims <- range(lapply(tracks,values))
  
  tracks_to_plot<-OverlayTrack(tracks,name=as.character(patient),
                               ylim=ylims)
  
  annotracks <- peak_tracks[combos]
  
  
  plotTracks(background.title="white",fontcolor.title="black",col.axis="black",
             showTitle=T, clip=FALSE,
             c(AnnotationTrack(peaks, name = NULL, fill="grey"),
               annotracks,
               tracks_to_plot),
             legend=legend,
             type="l",
             from=start,
             to = end,
             min.distance=10,
             chromosome=chr,
             lwd.baseline=1,
             col.baseline="black",
             windowSize =40,window=-1,
             main=paste(mcols(peaks)[mcols(peaks)$peak==peak,"id"],
                        mcols(peaks)[mcols(peaks)$peak==peak,"type"] ,
                        mcols(peaks)[mcols(peaks)$peak==peak,"event_type_chr"],"\n",
                        peak, 
                        sep=" "),
             ylim=ylims,cex.title=1,cex.axis=0.7, cex.main=1, add=TRUE,
             sizes=c(0.02,rep(0.02,length(combos)),1),
             detail="coverage")
}



col <- brewer.pal(name="Set1", n=9)
names(col) <- c("A","B","C","D","E","F","G","H","I")


plot_purity_C561 <- function(patient=patient){
  purity_plot <- coldata[coldata$patient==patient & !grepl("E",coldata$region),]
  purity_plot$region <- as.factor(as.character(purity_plot$region)) 
  ggplot(purity_plot, aes(region,purity, fill=region), ylim=c(0,1)) + 
    geom_boxplot() +
    theme_classic() +
    scale_y_continuous(limits = c(0,1.0), expand=c(0,0),labels = scales::percent) +
    theme(axis.line.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "none") +
    scale_fill_manual(values = col)
}


Pval_table_C561 <- function(peak=peak){
  NvA_logfc_F <- round(fc[peak,"C561.adenoma..F."], digits=2)
  NvA_pvalue_F <- p[peak,"C561.adenoma..F."]
  NvA_padj_F <- p_adj[peak,"C561.adenoma..F."]
  NvA_logfc_G <- round(fc[peak,"C561.adenoma..G."], digits=2)
  NvA_pvalue_G <- p[peak,"C561.adenoma..G."]
  NvA_padj_G <- p_adj[peak,"C561.adenoma..G."]
  subtissue_logfc_F <- subtissue_nf[subtissue_nf$peak==peak,"F_logfc"]
  subtissue_pvalue_F <- subtissue_nf[subtissue_nf$peak==peak,"F_pvalue"]
  subtissue_padj_F <- subtissue_nf[subtissue_nf$peak==peak,"F_padj"]
  subtissue_logfc_G <- subtissue_nf[subtissue_nf$peak==peak,"G_logfc"]
  subtissue_pvalue_G <- subtissue_nf[subtissue_nf$peak==peak,"G_pvalue"]
  subtissue_padj_G <- subtissue_nf[subtissue_nf$peak==peak,"G_padj"]
  results_plot <- rbind(NvA_pvalue,NvA_padj,subtissue_pvalue_F,subtissue_padj_F,subtissue_pvalue_G,subtissue_padj_G)
  results_plot<-as.data.frame(as.numeric(results_plot[,1]))
  sig<-ifelse(is.na(results_plot[,1]),"",
              ifelse(results_plot[,1] <= 0.0001 , " ****",
                     ifelse(results_plot[,1] <= 0.001 , " ***",
                            ifelse(results_plot[,1] <= 0.01 , " **",
                                   ifelse(results_plot[,1] <= 0.05 , " *", "")))))
  
  results_plot <- as.data.frame(sapply(sapply(results_plot, as.numeric), scientific))
  results_plot[,1] <- paste0(results_plot[,1],sig)
  results_plot<-rbind(results_plot,NvA_logfc,subtissue_logfc_F,subtissue_logfc_G)
  LogFC <- c(results_plot[7:9,1])
  P. <- c(results_plot[c(1,3,5),1])
  `Adjusted P.` <- c(results_plot[c(2,4,6),1])
  results_plot <- cbind(LogFC,P.,`Adjusted P.`)
  rownames(results_plot) <- c("Adenoma vs. normal","F vs.Cancer ~purity+type","G vs.Cancer ~purity+type")
  
  title <- textGrob(get_status_adenoma_C561(peak),
                    gp=gpar(fontsize=12))

  padding <- unit(5,"mm")
  t1<-tableGrob(results_plot, theme=ttheme_default(base_size=10))
  table <- gtable_add_rows(
    t1, 
    heights = grobHeight(title) + padding,
    pos = 0)
  table <- gtable_add_grob(
    table, 
    title, 
    1, 1, 1 ,ncol(table))
  grid.draw(table)
}


get_status_adenoma_C561 <- function(peak=peak){
  subtissue_padj_F <- subtissue_nf[subtissue_nf$peak==peak,"F_padj"]
  subtissue_padj_G <- subtissue_nf[subtissue_nf$peak==peak,"G_padj"]
  NvA_padj_F <- p_adj[peak,"C561.adenoma..F."]
  NvA_padj_G <- p_adj[peak,"C561.adenoma..G."]
  if(is.na(subtissue_padj_F)&is.na(subtissue_padj_G)&is.na(NvA_padj_F)&is.na(NvA_padj_G)){title <- "NA"
  } else if(!is.na(subtissue_padj_F)&subtissue_padj_F <0.05 | !is.na(subtissue_padj_G)& subtissue_padj_G <0.05 ){ title <- "Putatively Subclonal"
  } else if(sig_filter[peak,"C561.adenoma..F."]==TRUE | sig_filter[peak,"C561.adenoma..G."]==TRUE | sig_filter[peak,"C561.adenoma..H."]==TRUE){ title <- "Putatively Clonal"
  } else { title <- "NS"}
  return(title)
}

plot_info_C561 <- function(pdf_name, patient, peak, type){
  pdf(pdf_name, height=5, width=10.4)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(3, 2, widths = c(6,4.4))))
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  pushViewport(viewport(layout.pos.row=1:3,layout.pos.col=1))
  print(plot_horizon_C561(cov=cov,
                     patient=patient,
                     peak=peak,
                     smooth=30,
                     lwd=1.5,
                     legend=T,
                     type="l"), vp= vplayout(x=1:3,y=1))
  popViewport()
  pushViewport(viewport(layout.pos.row=1,layout.pos.col=2))
  print(Pval_table_C561(peak=peak),vp = vplayout(3, 2))
  popViewport()
  pushViewport(viewport(layout.pos.row=2:3,layout.pos.col=2))
  p1<-plot_purity_C561(patient=patient)
  print(p1,vp = vplayout(2:3,2))
  popViewport()
  graphics.off()
}





