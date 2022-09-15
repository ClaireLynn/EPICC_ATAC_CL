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
dir <- args[4] 

counts <- readRDS(args[1])
all_samp <- names(mcols(counts))
patient <- sapply(strsplit(basename(args[1]),"_"),`[`,1)
file <- sub(".rds","",basename(args[1]))

model <- as.formula(args[2])
reduced <- as.formula(args[3])

coldata <- readRDS("coldata.rds")
peaks <- readRDS("recurrent_peaks.rds")
peaks_filt <- peaks

# Remove samples with no purity info
coldata <- coldata[!(is.na(coldata$purity) | coldata$purity == 0), ]

# Remove samples with no purity info from counts
mcols(counts) <- mcols(counts)[,colnames(mcols(counts)) %in% rownames(coldata)]

# Select only coldata of this patient
coldata <- coldata[colnames(mcols(counts)),]
coldata$purity_short <- round(coldata$purity*100, digits=0)


# Make the DESeq2 objects
count_matrix <- as.matrix(mcols(counts))
rownames(count_matrix) <- paste0(seqnames(counts),":",start(counts),"-",end(counts))
regions <- as.character(unique(coldata$region))
count_matrix <- unique(count_matrix)
# Select only recurrent peaks
count_matrix_fil <- count_matrix[peaks$peak,]
count_matrix_fil <- unique(count_matrix_fil)

dds_fil <- DESeqDataSetFromMatrix(countData=count_matrix_fil,
         colData=coldata,
         design=model)


# Get normfacs from non-filtered peak set!
nfpeakSizeFactors <-readRDS(paste0("normfacs/",patient,"_nf_counts_patientSpecPeakSizeFactors_adenoma.rds"))


# Do DESEQ
dds_fil$sizeFactor <- nfpeakSizeFactors
dds_fil <- estimateDispersions(dds_fil)
dds_fil <- DESeq(dds_fil, test="LRT", reduced=reduced)
res <- results(dds_fil)

saveRDS(res,paste0(dir,"/",file,gsub(" ","",as.character(model)[2]),"_deseq_LRT.rds"))
