# This script performs DESEQ per patient, between regions
suppressMessages(library(DESeq2))

# Get arguments and coldata
args = commandArgs(trailingOnly=TRUE)
# directory name, depends on the test
dir <- args[4] 

patient <- args[1]
counts <- readRDS(paste0("counts_recurrent/",patient,"_recurrent_peaks_nf_counts.rds"))
all_samp <- names(mcols(counts))

file <- sub(".rds","",basename(counts))
# get model, one of: ~purity+sample_type+region ~sample_type+purity ~sample_type+region
model <- as.formula(args[2])
# get reduced model, one of: ~sample_type+purity ~sample_type 
reduced <- as.formula(args[3])

coldata <- readRDS("coldata.rds")
peaks <- readRDS("recurrent_peaks.rds")
peaks_filt <- peaks

# Remove samples with no purity info/ adenomas from coldata
coldata <- coldata[!(is.na(coldata$purity) | coldata$purity == 0 | coldata$subtissue == "adenoma"), ]

# Remove samples with no purity info/ adenomas from counts 
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

# Get normfacs from non-filtered peak set
nfpeakSizeFactors <-readRDS(paste0("normfacs/",patient,"_",args[5],"_counts_patientSpecPeakSizeFactors.rds"))

# Do DESEQ
dds_fil$sizeFactor <- nfpeakSizeFactors
dds_fil <- estimateDispersions(dds_fil)
dds_fil <- DESeq(dds_fil, test="LRT", reduced=reduced)
res <- results(dds_fil)

saveRDS(res,paste0(dir,"/",file,gsub(" ","",as.character(model)[2]),"_deseq_LRT.rds"))


