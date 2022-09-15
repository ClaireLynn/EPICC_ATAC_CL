# This script selects the peaks to test/plot between regions, based on the top 20 most recurrent:
# gained promoter, lost promoter, gained enhancer, lost enhancer
# Plus driver mutations with recurrence >=5
library(stringr)
library(tidyr)
library(dplyr)

# get genehancer information
tryCatch({
  genehancer_data = "genehancer_tracks" %>%
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

# load recurrence data
final_reccurence_count = readRDS("final_reccurence_data.rds")

# get recurrence count in cancer
final_reccurence_count$sig_matrix[is.na(final_reccurence_count$sig_matrix)]<- FALSE
s_matrix = final_reccurence_count$sig_matrix
s_matrix[is.na(s_matrix)] = FALSE
cancer_s_matrix <-s_matrix[,grepl("pure",colnames(s_matrix))]
adenoma_s_matrix <-s_matrix[,grepl("adenoma",colnames(s_matrix))]
final_reccurence_count$summary$recurrence_cancer <- rowSums(cancer_s_matrix)
final_reccurence_count$summary$recurrence_adenoma <- rowSums(adenoma_s_matrix)

peaks_to_plot = 
  final_reccurence_count$summary %>% 
  mutate(peak=as.character(peak)) %>%
  arrange(-recurrence_cancer,id
  ) %>% 
  # exclude non consistent changes and keep top 5 gain/losses for enhancers/promotors
  filter(!is.na(event_type)) %>%
  filter(event_type!="gain/loss")  %>%
  # split by event type and select most reccurent ones
  split(., list(.$type, .$event_type)) %>% 
  lapply(head, n=20) %>% 
  do.call(what=rbind) %>% 
  # sort by reccurence and add addtional data to label
  mutate(id=as.character(id)) %>% 
  mutate(id=ifelse(type == "enhancer", longer_enh_label[id], id)) %>% 
  mutate(id=factor(peak, peak[!duplicated(peak)], make.unique(id[!duplicated(peak)])))

# get driver list
drivers <-readLines("IntOGen-DriverGenes_COREAD.txt")
drivers <- drivers[order(drivers)]

summary <- final_reccurence_count$summary
summary$symbol <- summary$id
summary=
  summary %>%
  mutate(symbol=ifelse(type == "enhancer", longer_enh_label[symbol], symbol)) %>%
  filter(!is.na(event_type)) %>%
  filter(!(event_type=="gain/loss")) %>%
  filter(!(recurrence<5)) %>%
  filter(!(peak %in% peaks_to_plot$peak)) 

# get the drivers
Pattern = paste0("(^| |,)",drivers,"($| |,)", collapse="|")
driver_res =summary[grepl(Pattern, summary$symbol),]
driver_res=driver_res  %>% group_by(id) %>%  dplyr::slice(which.max(recurrence))
driver_res =driver_res[order(driver_res$type,driver_res$event_type,-driver_res$recurrence),]

extra_peaks_to_plot = 
  final_reccurence_count$summary %>% 
  mutate(peak=as.character(peak)) %>% 
  filter(peak %in% driver_res$peak) %>%
  mutate(id=as.character(id)) %>% 
  mutate(id=ifelse(type == "enhancer", longer_enh_label[id], id)) %>% 
  mutate(id=factor(peak, peak[!duplicated(peak)], make.unique(id[!duplicated(peak)])))

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
saveRDS(peaks,"recurrent_peaks.rds")
