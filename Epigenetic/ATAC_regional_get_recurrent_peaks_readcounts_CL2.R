# This script gets the per-sample read counts for the recurrent peaks

library(GenomicRanges)
library(stringr)
if (!require(optparse)) stop("Package 'optparse' missing\n.")

source("~/projects/EPICC/ATAC/4-excluded_samples.R")
excluded_samples <- str_match(excluded_samples, "EPICC_\\s*(.*?)\\s*_C1$")[,2]
excluded_samples <- c(excluded_samples, "C537_D1_G8")
option_list = list(
        make_option(c("-p","--patient"), type="character", default=NULL, metavar="character"),
        make_option(c("-r","--regions"), type="character", default=NULL, metavar="character"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

regions = opt$regions
peaks <- readRDS(regions)
peak_file <- sub(".rds","",basename(regions))

patient = opt$patient

nf_cutsites <- list.files(path = "~/projects/EPICC/ATAC/cutsites",
                          pattern = ".*[1-9]_nucleosome_free-GRCh38-proper_pairs-shifted.bed.gz$",
			  full.names = TRUE,
			  recursive=TRUE)

nf_cutsites <- nf_cutsites[grepl(patient,nf_cutsites)]

samples <- str_match(nf_cutsites, "EPICC_\\s*(.*?)\\s*_C1_nucleosome_free-GRCh38-proper_pairs-shifted.bed.gz")[,2]
exclude <- samples %in% excluded_samples
samples <- samples[!exclude]

nf_cutsites <- nf_cutsites[!exclude]

gr <- lapply(1:length(samples), function(s) {
df<- read.table(nf_cutsites[s])
colnames(df) <- c("chr","start","end","strand")
gr <- makeGRangesFromDataFrame(df,
                         keep.extra.columns=FALSE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         starts.in.df.are.0based=FALSE)
gr <- resize(gr, width = 100, fix = "center")

peak_counts <- countOverlaps(peaks,gr)
print(paste0(samples[s]," peak counts done"))

return(peak_counts)
})

values(peaks) <- gr
rm(gr)
colnames(values(peaks)) <- samples
saveRDS(peaks,file=paste(patient,peak_file,"nf_counts.rds",sep="_"))

