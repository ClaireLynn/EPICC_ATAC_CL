# This script gets the normalisation factors for DESEQ (using counts in all peaks)

suppressMessages(library(DESeq2))

# Write sessionInfo
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

# Get arguments and coldata
args = commandArgs(trailingOnly=TRUE)
dir <- "normfacs"

patient <- args[1]
# get counts, just for sample info
counts <- readRDS(paste0("counts_recurrent/",patient,"_recurrent_peaks_nf_counts.rds"))
all_samp <- names(mcols(counts))
file <- paste0(patient,"_nf_counts_patientSpecPeakSizeFactors.rds")

model <- as.formula("~purity+sample_type")

coldata <- readRDS("coldata.rds")

# Remove samples with no purity info/ adenomas from coldata
coldata <- coldata[!(is.na(coldata$purity) | coldata$purity == 0 | coldata$subtissue == "adenoma") , ]

# Remove samples with no purity info/ adenomas from counts 
mcols(counts) <- mcols(counts)[,colnames(mcols(counts)) %in% rownames(coldata)]

# Select only coldata of this patient
coldata <- coldata[colnames(mcols(counts)),]
coldata$purity_short <- round(coldata$purity*100, digits=0)

# Get dataset for normfacs
normdata <- readRDS("peak_coverage_data_nucleosome_free_filtered.rds")
count_matrix <- normdata$data
count_matrix<-count_matrix[grepl(paste0(patient,".region_A|",
                                        patient,".region_B|",
                                        patient,".region_C|",
                                        patient,".region_D|",
                                        patient,".region_F|",
                                        patient,".region_G|",
                                        patient,".region_H|",
                                        patient,".region_I|"),
                                 normdata$peak_calls$name),]

a <- unlist(lapply(strsplit(colnames(count_matrix),"-"),`[[`,2))
b <- unlist(lapply(strsplit(a,"_"),`[[`,2))
c <- unlist(lapply(strsplit(a,"_"),`[[`,3))
d <- unlist(lapply(strsplit(a,"_"),`[[`,4))

colnames(count_matrix) <- paste(b,c,d,sep="_")
count_matrix <- count_matrix[,colnames(count_matrix) %in% rownames(coldata)]

dds <- DESeqDataSetFromMatrix(countData=count_matrix,
         colData=coldata,
         design=model)

# Get normfacs from  peak set
peakSizeFactors <-estimateSizeFactors(dds)$sizeFactor

saveRDS(peakSizeFactors,paste0(dir,"/",file))
