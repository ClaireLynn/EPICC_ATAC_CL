suppressMessages(library(DESeq2))
suppressMessages(library(edgeR))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(ggpubr))
suppressMessages(library(RColorBrewer))
library(csaw)
set.seed(123)

# Write sessionInfo
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

# Get arguments and coldata
args = commandArgs(trailingOnly=TRUE)
dir <- args[7] 

counts <- readRDS(args[1])
all_samp <- names(mcols(counts))
patient <- sapply(strsplit(basename(args[1]),"_"),`[`,1)
file <- paste0(patient,"_nf_counts_patientSpecPeakSizeFactors_adenoma.rds")
model <- as.formula(args[2])
reduced <- as.formula(args[3])

coldata <- readRDS("coldata.rds")

# Remove samples with no purity info/ adenomas from coldata
coldata <- coldata[!(is.na(coldata$purity) | coldata$purity == 0) , ]

# Remove samples with no purity info/ adenomas from counts 
mcols(counts) <- mcols(counts)[,colnames(mcols(counts)) %in% rownames(coldata)]

# Select only coldata of this patient
coldata <- coldata[colnames(mcols(counts)),]
coldata$purity_short <- round(coldata$purity*100, digits=0)


# Get dataset for normfacs
normdata <- readRDS(args[4])
count_matrix <- normdata$data
count_matrix<-count_matrix[grepl(paste0(patient,".region_A|",patient,".region_B|",patient,".region_C|",patient,".region_D|",patient,".region_F|",patient,".region_G|",patient,".region_H|",patient,".region_I|"),normdata$peak_calls$name),]

a <- unlist(lapply(strsplit(colnames(count_matrix),"-"),`[[`,2))
b <- unlist(lapply(strsplit(a,"_"),`[[`,2))
c <- unlist(lapply(strsplit(a,"_"),`[[`,3))
d <- unlist(lapply(strsplit(a,"_"),`[[`,4))

colnames(count_matrix) <- paste(b,c,d,sep="_")
count_matrix <- count_matrix[,colnames(count_matrix) %in% rownames(coldata)]

dds <- DESeqDataSetFromMatrix(countData=count_matrix,
         colData=coldata,
         design=model)

# Get normfacs from  peak set!
peakSizeFactors <-estimateSizeFactors(dds)$sizeFactor

saveRDS(peakSizeFactors,paste0(dir,"/",file))

