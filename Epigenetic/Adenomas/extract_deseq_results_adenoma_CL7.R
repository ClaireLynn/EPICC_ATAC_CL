library(stringr)
peaks <- readRDS("~/EPICC_data/DESEQ_final/recurrent_peaks.rds")
cov <- readRDS("~/EPICC_data/recurrent_peaks_final/cpm_peaknorm_cov_nf_1bp.rds")

clonal <-"Putatively Clonal"
subclonal <-"Putatively Subclonal"

figure <-c("C561", "C532", "C551", "C552", "C531", "C538", "C539", "C537", "C525", "C528",
           "C524", "C530", "C550", "C544", "C554", "C518", "C562", "C548","C536", "C559",
           "C516", "C543", "C542", "C555")

coldata <- readRDS("~/EPICC_data/DESEQ_final/coldata.rds")
coldata_region <- readRDS("~/EPICC_data/DESEQ_final/coldata_region.rds")
combos <- rownames(coldata_region)
patients <- unique(coldata_region$patient)

patients <- unique(coldata_region[coldata_region$subtissue == "adenoma","patient"])

patients <- patients[patients %in% figure]

setwd("~/EPICC_data/DESEQ_final/")
files<-list.files(path ="effect_of_subtissue_nfreads_LRT",pattern="LRT.rds", full.names = TRUE)
files <- files[grepl(paste(figure,collapse="|"),files)]
subtissue_nf <- lapply(1:length(files), function(i){
  pvalue <- as.vector(readRDS(files[i])$pvalue)
  padj <- as.vector(readRDS(files[i])$padj)
  logfc <- -round(as.vector(readRDS(files[i])$log2FoldChange),digits=2)
  pvalue <- as.data.frame(cbind(row.names(readRDS(files[i])),pvalue,padj,logfc))
})
for(i in 1:length(subtissue_nf)){
  colnames(subtissue_nf[[i]]) <- c("peak",
                                   paste0(unlist(str_split(unlist(str_split(files[i],"/"))[2], "_"))[1],"_pvalue"),
                                   paste0(unlist(str_split(unlist(str_split(files[i],"/"))[2], "_"))[1],"_padj"),
                                   paste0(unlist(str_split(unlist(str_split(files[i],"/"))[2], "_"))[1],"_logfc"))
}
subtissue_nf<-Reduce(function(x, y) merge(x, y, by="peak", all=T), subtissue_nf)



# Get the bulk test significance!

# Get the bulk test significance!
sig_filter<- readRDS("~/git_projects/EPICC_ATAC/final_reccurence_data.rds")$sig_matrix[peaks$peak,]
fc <- readRDS("~/git_projects/EPICC_ATAC/final_reccurence_data.rds")$fc[peaks$peak,]
p <- readRDS("~/git_projects/EPICC_ATAC/final_reccurence_data.rds")$p[peaks$peak,]
p_adj <- readRDS("~/git_projects/EPICC_ATAC/final_reccurence_data.rds")$p_adj[peaks$peak,]

