library(ggplot2)
library(cowplot)
library(tidyr)
library(dplyr)

tryCatch({
  genehancer_data = 
    "~/git_projects/EPICC_ATAC/genehancer_tracks" %>%
    list.files("[.]bed[.]gz", full.names=TRUE) %>% 
    (function(x) x[grepl("DoubleElite", x)]) %>%
    magrittr::set_names(., gsub("(DoubleElite)?[.]bed[.]gz","", basename(.))) %>%
    lapply(readr::read_tsv)
  longer_enh_label = 
    with(genehancer_data$geneHancerInteractions, {
      e2g = sapply(split(geneName, geneHancerIdentifier), paste, collapse=", ")
      long_label = paste0(e2g, " (", names(e2g),")")
      names(long_label) = names(e2g)
      return(long_label)
    })
}, error=function(e) print("Couldn't load genehancer dataset."))

longer_enh_label<-unlist(lapply(lapply(strsplit(longer_enh_label,","),tail ,n=3),paste,collapse="," ))

tryCatch({
  genehancer_data = 
    "~/git_projects/EPICC_ATAC/genehancer_tracks" %>%
    list.files("[.]bed[.]gz", full.names=TRUE) %>% 
    (function(x) x[grepl("DoubleElite", x)]) %>%
    magrittr::set_names(., gsub("(DoubleElite)?[.]bed[.]gz","", basename(.))) %>%
    lapply(readr::read_tsv)
  enh_symbols = 
    with(genehancer_data$geneHancerInteractions, {
      e2g = sapply(split(geneName, geneHancerIdentifier), paste, collapse=", ")
      long_label = paste0(e2g)
      names(long_label) = names(e2g)
      return(long_label)
    })
}, error=function(e) print("Couldn't load genehancer dataset."))

theme_set(theme_cowplot())
subclonal_matrix<-readRDS("~/git_projects/EPICC_ATAC/subclonal_matrix.rds")
final_reccurence_count = readRDS("~/git_projects/EPICC_ATAC/final_reccurence_data.rds")
rnaseq_test_results = readRDS("~/git_projects/rnaseq_test_results_all_final.rds")
rnaseq_test_results$p_adj_deseq <- p.adjust(rnaseq_test_results$p_deseq, "fdr")
fig_dir = "~/subclonal_redo/"
# significant fold change
figure <-c("C561", "C532", "C551", "C552", "C531", "C538", "C539", "C537", "C525", "C528",
           "C524", "C530", "C550", "C544", "C554", "C518", "C562", "C548","C536", "C559",
           "C516", "C543", "C542", "C525")



final_reccurence_count$p_adj<- final_reccurence_count$p_adj[,!grepl("C547", colnames(final_reccurence_count$p_adj))]
final_reccurence_count$fc<- final_reccurence_count$fc[,!grepl("C547", colnames(final_reccurence_count$fc))]
final_reccurence_count$sig_matrix<- final_reccurence_count$sig_matrix[,!grepl("C547", colnames(final_reccurence_count$sig_matrix))]

#final_reccurence_count$p_adj<- final_reccurence_count$p_adj[,!grepl("C561.adenoma..H.", colnames(final_reccurence_count$p_adj))]
#final_reccurence_count$fc<- final_reccurence_count$fc[,!grepl("C561.adenoma..H.", colnames(final_reccurence_count$fc))]
#final_reccurence_count$sig_matrix<- final_reccurence_count$sig_matrix[,!grepl("C561.adenoma..H.", colnames(final_reccurence_count$sig_matrix))]

final_reccurence_count$sig_matrix[is.na(final_reccurence_count$sig_matrix)]<- FALSE
s_matrix = final_reccurence_count$sig_matrix
s_matrix[is.na(s_matrix)] = FALSE
cancer_s_matrix <-s_matrix[,grepl("pure",colnames(s_matrix))]
adenoma_s_matrix <-s_matrix[,grepl("adenoma",colnames(s_matrix))]
final_reccurence_count$summary$recurrence_cancer <- rowSums(cancer_s_matrix)
final_reccurence_count$summary$recurrence_adenoma <- rowSums(adenoma_s_matrix)

fc_sign = sign(final_reccurence_count$fc[rownames(s_matrix),colnames(s_matrix)])
fc_sign[is.na(s_matrix) | !s_matrix] = 0
# correct direction of change
n_losses = apply(fc_sign < 0, 1, sum)
n_gains = apply(fc_sign > 0, 1, sum)
change_type = case_when(n_losses < n_gains~1, n_losses > n_gains ~ -1, TRUE ~ 0)
names(change_type) = names(n_losses)
for (i in seq_len(NCOL(fc_sign))) { 
  fc_sign[which(fc_sign[,i] != change_type),i] = 0  
}

min_p_rnaseq = with(rnaseq_test_results, tapply(p_deseq, peak.peak, min))
min_p_rnaseq = with(rnaseq_test_results, tapply(p_deseq, peak.peak, min))
rnaseq_test_results <-rnaseq_test_results[!is.na(rnaseq_test_results$p_adj_deseq),]
min_p_rnaseq_fdr = with(rnaseq_test_results, tapply(p_adj_deseq, peak.peak, min))

# 
peaks_to_plot = 
  final_reccurence_count$summary %>% 
  mutate(peak=as.character(peak)) %>% 
  mutate(event_type = change_type[peak]) %>% 
  mutate(event_type_chr = factor(event_type, c(-1,1), c("loss","gain"))) %>% 
  mutate(p_rna_correlation = min_p_rnaseq_fdr[peak]) %>% 
 #mutate(used_before=peak %in% rownames(sc_cancer)) %>% 
 arrange(-recurrence_cancer,id
         #-used_before
         ) %>% 
  # exclude non consistent changes and keep top 5 gain/losses for enhancers/promotors
  filter(!is.na(event_type_chr)) %>%
  # split by event type and select most reccurent ones
  split(., list(.$type, .$event_type_chr)) %>% 
  lapply(head, n=20) %>% 
  do.call(what=rbind) %>% 
  # sort by reccurence and add addtional data to label
  #arrange(-recurrence)  %>% 
  mutate(id=as.character(id)) %>% 
  mutate(id=ifelse(type == "enhancer", longer_enh_label[id], id)) %>% 
  mutate(id=factor(peak, peak[!duplicated(peak)], make.unique(id[!duplicated(peak)])))

extra_peaks_to_plot = 
  final_reccurence_count$summary %>% 
  mutate(peak=as.character(peak)) %>% 
  mutate(event_type = change_type[peak]) %>% 
  mutate(event_type_chr = factor(event_type, c(-1,1), c("loss","gain"))) %>% 
  mutate(p_rna_correlation = min_p_rnaseq_fdr[peak]) %>% 
  #mutate(used_before=peak %in% rownames(sc_cancer)) %>% 
  arrange(-recurrence_cancer,id
          #-used_before
  ) %>% 
  # exclude non consistent changes and keep top 5 gain/losses for enhancers/promotors
  filter(!is.na(event_type_chr)) %>%
  mutate(id=as.character(id)) %>% 
  mutate(id=ifelse(type == "enhancer", longer_enh_label[id], id)) %>% 
  mutate(id=factor(peak, peak[!duplicated(peak)], make.unique(id[!duplicated(peak)])))


# merged <- merge(final_reccurence_count$sig_matrix,final_reccurence_count$fc, by=0, suffixes=c("_sig","_fc"))
# rownames(merged) <- merged$Row.names
# merged$Row.names <- NULL
# merged <- merged[,order(colnames(merged))]
# merged <- merge(extra_peaks_to_plot,merged, by=0)
# rownames(merged) <- merged$Row.names
# merged$Row.names <- NULL
# merged <- merged[,c(2,6,3,4,10,5,8,9,11,12:ncol(merged))]
# 
# frac_up <- vector()
# for(pat in patients){
# sig_name<-paste0(pat,".pure_sig")
# fc_name<-paste0(pat,".pure_fc")
# mergedd <- merged[merged[,sig_name]==TRUE,]
# frac_up[pat] <-sum(mergedd[mergedd[,sig_name],fc_name]>0)/nrow(mergedd)
# }
# 
# totals <- vector()
# for(pat in patients){
#   sig_name<-paste0(pat,".pure_sig")
#   fc_name<-paste0(pat,".pure_fc")
#   mergedd <- merged[merged[,sig_name]==TRUE,]
#   totals[pat] <- nrow(mergedd)
# }
# 
# 
# frac_up2 <- vector()
# merged2 <- merged[merged$recurrence_cancer>4,]
# for(pat in patients){
#   sig_name<-paste0(pat,".pure_sig")
#   fc_name<-paste0(pat,".pure_fc")
#   mergedd <- merged2[merged2[,sig_name]==TRUE,]
#   frac_up2[pat] <-sum(mergedd[mergedd[,sig_name],fc_name]>0)/nrow(mergedd)
# }
# 
# write.csv(merged,"~/subclonal_redo/merged_results.csv")


#extra_peaks_from_andrea <- unique(c("chr5:160213457-160213957",
#                             "chr20:3168052-3168552",
#                             "chr11:62542695-62543195",
#                             "chr7:101237662-101238162",
#                             "chr20:31608508-31609008",
#                             "chr1:110390560-110391060",
#                             "chr1:151992421-151992921",
#                             "chr1:209790724-209791224",
#                             "chr3:119781791-119782291",
#                             "chr3:126179733-126180233",
#                             "chr5:141223078-141223578", driver_res$peak[-c(2,6)]))

extra_peaks_from_andrea <- driver_res$peak[-c(2,6)]

extra_peaks_to_plot <- extra_peaks_to_plot[extra_peaks_to_plot$peak %in% extra_peaks_from_andrea,]
#extra_peaks_to_plot$id <- gsub("\\.1","",extra_peaks_to_plot$id )
peaks_to_plot <- rbind(peaks_to_plot,extra_peaks_to_plot)

#making the peaks for subclonal analysis
chr <- unlist(lapply(str_split(peaks_to_plot$peak, ":"), `[[`, 1))
start <- unlist(lapply(str_split(lapply(str_split(peaks_to_plot$peak, ":"), `[[`, 2), "-"),  `[[`, 1))
end <- unlist(lapply(str_split(lapply(str_split(peaks_to_plot$peak, ":"), `[[`, 2), "-"),  `[[`, 2))
df <- as.data.frame(cbind(chr,start,end,peaks_to_plot))
row.names(df) <- NULL
peaks <- makeGRangesFromDataFrame(df,
                                  keep.extra.columns=TRUE,
                                  ignore.strand=FALSE,
                                  seqinfo=NULL,
                                  starts.in.df.are.0based=FALSE)
saveRDS(peaks,"~/git_projects/EPICC_ATAC/recurrent_peaks.rds")


#
effect_group_labels = # final label for plot
  c("Gained Promoters"="Gained Promoters", 
    "Lossed Promoters"="Lost Promoters", 
    "Gained Enhancers"="Gained Enhancers", 
    "Lossed Enhancers"="Lost Enhancers")
plot_data = # actuall statistics for the peaks 
  fc_sign[peaks_to_plot$peak,] %>% 
  mutate(peak=rownames(.)) %>% 
  reshape2::melt(id.vars="peak") %>% 
  magrittr::set_colnames(c("peak","patient","effect")) %>% 
  merge(peaks_to_plot, all.x=TRUE, by="peak") %>% 
  mutate(effect_group=paste0(Hmisc::capitalize(as.character(event_type_chr)), "ed ", Hmisc::capitalize(type), "s")) %>% 
  mutate(effect_group=factor(effect_group, names(effect_group_labels), effect_group_labels)) %>% 
  mutate(status=ifelse(grepl("C516|C536.pure|C548.pure|C552.pure|C518.pure", patient), "MSI", "MSS")) %>%
  mutate(status=factor(status, c("MSS","MSI"), labels = c("MSS","MSI"))) %>%
  mutate(tissue=ifelse(grepl("[.]pure", patient), "cancer", "adenoma")) %>%
  mutate(patient=gsub("[.]adenoma", "", gsub("[.]pure", "", patient))) %>% 
  mutate(patient=gsub("[.]\\.", " (", patient)) %>% 
  mutate(patient=gsub("\\.", ")", patient))   %>% 
  mutate(effect = ifelse(effect == c("gain"=1, "loss"=-1)[as.character(event_type_chr)], 1, NA)) %>% 
  mutate(tissue=factor(tissue, c("adenoma","cancer"), labels = c("Adenoma","Cancer"))) 

# set order of peaks and patients 
order_ids = rev(peaks_to_plot$id)
wh = !is.na(plot_data$effect) & plot_data$tissue == "Cancer"
check_order <- table(plot_data$patient[wh])
check_order <- check_order[order(check_order, decreasing=T)]
#ord_pats = tapply(2^match(plot_data$id[wh], order_ids), plot_data$patient[wh], sum)
#ord_pats[is.na(ord_pats)] = 0
#ord_pats = unique(c(names(sort(ord_pats, decreasing = TRUE)), plot_data$patient))
ord_pats = c(names(check_order),"C561 (F)","C561 (G)","C561 (H)"
             )
plot_data = 
  plot_data %>% 
  mutate(id=factor(id, order_ids, ordered = TRUE)) %>%
  mutate(patient=factor(patient, ord_pats, ordered = TRUE))
# RNA-seq tests to add to plot:
plot_data_rna = 
  dplyr::filter(plot_data, p_rna_correlation < 0.01 & !is.na(effect))

# change ids to include symbol for enhancers


#plot_data[plot_data$type=="enhancer","new_id"] <- paste0(plot_data[plot_data$type=="enhancer","symbol"],
 #                                                        " ","(",plot_data[plot_data$type=="enhancer","id"],")")

# create the plot
plot_heatmap =
  plot_data %>%
  mutate(effect=ifelse(is.na(effect), 0, effect)) %>%
  ggplot(aes(x=patient, y=id, fill=event_type_chr, alpha=effect*0.9)) + 
  # geom layers:
  geom_tile(color="gray0", size=0.2) + 
  geom_point(data=plot_data_rna, aes(color="Association with RNA-seq"), size=0.8) + 
  facet_grid(effect_group~status+tissue, labeller=labeller(.default=Hmisc::capitalize), scales="free", space="free") + 
  # scales and legends
  scale_color_manual(values=c("Association with RNA-seq"="gray10")) + 
  scale_alpha_identity() + 
  scale_fill_brewer(na.value="white", palette = "Set1", direction = -1) + 
  guides(color = guide_legend(title.position="top", title.hjust = 0.5), fill = FALSE) +
  # modify theme:
  theme_cowplot() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + 
  theme(legend.position="bottom", legend.box="horizontal") +
  theme(strip.text = element_text(size=11), strip.background = element_blank()) + 
  # labels
  xlab("") + ylab("") + labs(color="", fill="Fraction of samples mutated") +
  theme(strip.background=element_rect(fill=NA))
ofile = file.path(fig_dir, "heatmap_recurrent_epigenetic_changes_plusdrivers_RNA_deseqpadj0.01_final.pdf")
ggsave(ofile, plot_heatmap,width=10.5, height=17.1)




## second plot with external subclonal annotations
#sc_data = list(Cancer=subclonal_matrix, Adenoma=sc_adenoma)
sc_data = list(Cancer=subclonal_matrix)
plot_heatmap_sc = plot_heatmap
plot_heatmap_sc$data$effect = 
  sc_status = sapply(seq_len(NROW(plot_heatmap$data)), function(i) {
    idx_peak = as.character(plot_heatmap$data$peak[i])
    idx_pat = as.character(plot_heatmap$data$patient[i])
    idx_tissue = as.character(plot_heatmap$data$tissue[i])
    if (idx_pat %in% colnames(sc_data[[idx_tissue]]) & 
        idx_peak %in% rownames(sc_data[[idx_tissue]])) 
    {
      is_sc = sc_data[[idx_tissue]][idx_peak, idx_pat]
    } else {
      is_sc = FALSE
    }
    if (is_sc) return(0.4)
    return(plot_heatmap$data$effect[i])
  })
ofile = file.path(fig_dir, "heatmap_recurrent_epigenetic_changes_plus_sc_drivers_RNA_deseqpadj0.01_final.pdf")
ggsave(ofile, plot_heatmap_sc, width=10.5, height=17.1)

