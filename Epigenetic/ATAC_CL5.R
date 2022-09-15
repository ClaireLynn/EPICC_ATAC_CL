
library(stringr)
peaks_original <- readRDS("recurrent_peaks.rds")
# resize peaks to pad the track
peaks <- resize(peaks_original, width = 2000, fix = "center")
peaks_df <- mcols(peaks)

clonal <-"Putatively Clonal"
subclonal <-"Putatively Subclonal"

figure <-c("C561", "C532", "C551", "C552", "C531", "C538", "C539", "C537", "C525", "C528",
           "C524", "C530", "C550", "C544", "C554", "C518", "C562", "C548","C536", "C559",
           "C516", "C543", "C542", "C555")

coldata <- readRDS("coldata.rds")
coldata_region <- readRDS("coldata_region.rds")
combos <- rownames(coldata_region)
patients <- unique(coldata_region$patient)
patients <- figure

# Get the bulk test significance
sig_filter<- readRDS("final_reccurence_data.rds")$sig_matrix[peaks$peak,]
fc <- readRDS("final_reccurence_data.rds")$fc[peaks$peak,]
p <- readRDS("final_reccurence_data.rds")$p[peaks$peak,]
p_adj <- readRDS("final_reccurence_data.rds")$p_adj[peaks$peak,]

# Get DESEQ results
files<-list.files(path ="effect_of_purity_only_nfreads_LRT",pattern="LRT.rds", full.names = TRUE)
files <- files[grepl(paste(figure,collapse="|"),files)]
purity_only_nf <- lapply(1:length(files), function(i){
  pvalue <- as.vector(readRDS(files[i])$pvalue)
  padj <- as.vector(readRDS(files[i])$padj)
  pvalue <- as.data.frame(cbind(row.names(readRDS(files[i])),pvalue,padj))
})
for(i in 1:length(purity_only_nf)){
  colnames(purity_only_nf[[i]]) <- c("peak",
                                     paste0(unlist(str_split(unlist(str_split(files[i],"/"))[2], "_"))[1],"_purity_only_nf_pvalue"),
                                     paste0(unlist(str_split(unlist(str_split(files[i],"/"))[2], "_"))[1],"_purity_only_nf_padj"))
}
purity_only_nf<-Reduce(function(x, y) merge(x, y, by="peak", all=T), purity_only_nf)


files<-list.files(path ="effect_of_region_nfreads_LRT",pattern="LRT.rds", full.names = TRUE)
files <- files[grepl(paste(figure,collapse="|"),files)]
region_nf <- lapply(1:length(files), function(i){
  pvalue <- as.vector(readRDS(files[i])$pvalue)
  padj <- as.vector(readRDS(files[i])$padj)
  pvalue <- as.data.frame(cbind(row.names(readRDS(files[i])),pvalue,padj))
})
for(i in 1:length(region_nf)){
  colnames(region_nf[[i]]) <- c("peak",
                                paste0(unlist(str_split(unlist(str_split(files[i],"/"))[2], "_"))[1],"_region_nf_pvalue"),
                                paste0(unlist(str_split(unlist(str_split(files[i],"/"))[2], "_"))[1],"_region_nf_padj"))
}
region_nf<-Reduce(function(x, y) merge(x, y, by="peak", all=T), region_nf)


files<-list.files(path ="effect_of_region_only_nfreads_LRT",pattern="LRT.rds", full.names = TRUE)
files <- files[grepl(paste(figure,collapse="|"),files)]
region_only_nf <- lapply(1:length(files), function(i){
  pvalue <- as.vector(readRDS(files[i])$pvalue)
  padj <- as.vector(readRDS(files[i])$padj)
  pvalue <- as.data.frame(cbind(row.names(readRDS(files[i])),pvalue,padj))
})
for(i in 1:length(region_only_nf)){
  colnames(region_only_nf[[i]]) <- c("peak",
                                     paste0(unlist(str_split(unlist(str_split(files[i],"/"))[2], "_"))[1],"_region_only_nf_pvalue"),
                                     paste0(unlist(str_split(unlist(str_split(files[i],"/"))[2], "_"))[1],"_region_only_nf_padj"))
}
region_only_nf<-Reduce(function(x, y) merge(x, y, by="peak", all=T), region_only_nf)

results<-Reduce(function(x, y) merge(x, y, by="peak", all=T), list(mcols(peaks),
                                                                   region_only_nf,
                                                                   purity_only_nf,
                                                                   region_nf))

col <- brewer.pal(name="Set1", n=9)
names(col) <- c("A","B","C","D","E","F","G","H","I")


plot_purity <- function(patient=patient,peak=peak){
purity_plot <- coldata[coldata$patient==patient & grepl("A|B|C|D|C542_F",rownames(coldata)),]
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
results_plot <- results[results$symbol==peak, grepl(patient,colnames(results))]
CvN_res <- CvN[[patient]][CvN[[patient]]$symbol==peak,]
results_plot$CvN_adj <- CvN_res$PValue 
results_plot<-as.data.frame(sapply(results_plot, as.numeric))
sig<-ifelse(results_plot[,1] <= 0.0001 , " ****",
            ifelse(results_plot[,1] <= 0.001 , " ***",
                   ifelse(results_plot[,1] <= 0.01 , " **",
                          ifelse(results_plot[,1] <= 0.05 , " *", ""))))

results_plot <- as.data.frame(sapply(sapply(results_plot, as.numeric), scientific))

rownames(results_plot) <- c("~region","~purity","~purity+region","Normal Vs. Cancer")
colnames(results_plot) <- "P.Value"

results_plot$`P.Value` <- paste0(results_plot$`P.Value`,sig)
title <- textGrob(paste0(ifelse(sig[3]=="","Putatively Clonal","Putatively Subclonal"),
                         " (logFC = ",
                         round(CvN_res$logFC,digits = 2),
                         ")"),
                  gp=gpar(fontsize=14))
padding <- unit(5,"mm")
t1<-tableGrob(results_plot)
table <- gtable_add_rows(
  t1, 
  heights = grobHeight(title) + padding,
  pos = 0)
table <- gtable_add_grob(
  table, 
  title, 
  1, 1, 1, ncol(table))
grid.draw(table)
}

get_status <- function(patient=patient,peak=peak){
  results_plot <- results[results$symbol==peak, grepl(patient,colnames(results))]
  results_plot<-as.data.frame(sapply(results_plot, as.numeric))
  results_plot[is.na(results_plot)] <- 1
  sig<-ifelse(results_plot[,1] <= 0.0001 , "****",
              ifelse(results_plot[,1] <= 0.001 , "***",
                     ifelse(results_plot[,1] <= 0.01 , "**",
                            ifelse(results_plot[,1] <= 0.05 , "*", ""))))
  title <- ifelse(sig[3]=="","Clonal","Subclonal")
  return(title)
}

