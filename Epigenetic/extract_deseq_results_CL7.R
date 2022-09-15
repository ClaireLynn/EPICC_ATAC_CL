library(stringr)
setwd("/Users/clynn01/EPICC_data/DESEQ_final")
peaks_original <- readRDS("~/git_projects/EPICC_ATAC/recurrent_peaks.rds")
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

# Get the bulk test significance!
sig_filter<- readRDS("~/git_projects/EPICC_ATAC/final_reccurence_data.rds")$sig_matrix[peaks$peak,]
fc <- readRDS("~/git_projects/EPICC_ATAC/final_reccurence_data.rds")$fc[peaks$peak,]
p <- readRDS("~/git_projects/EPICC_ATAC/final_reccurence_data.rds")$p[peaks$peak,]
p_adj <- readRDS("~/git_projects/EPICC_ATAC/final_reccurence_data.rds")$p_adj[peaks$peak,]

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

saveRDS(results,"all_deseq_results.rds")