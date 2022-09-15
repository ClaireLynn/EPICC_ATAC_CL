# This script gets the coverage for recurrent peaks per region per patient
library(GenomicFeatures)
library(GenomicRanges)

coldata <- readRDS("coldata.rds")
coldata_region <- readRDS("coldata_region.rds")

figure <-c("C561", "C532", "C551", "C552", "C531", "C538", "C539", "C537", "C525", "C528",
           "C524", "C530", "C550", "C544", "C554", "C518", "C562", "C548","C536", "C559",
           "C516", "C543", "C542", "C555")

coldata_region <-  coldata_region[coldata_region$patient %in% figure,]
coldata <-  coldata[coldata$patient %in% figure,]

#get recurrent peaks as GRanges
peaks <- readRDS("recurrent_peaks.rds")
resizePeaks <- resize(peaks, width = 2000, fix = "center")
bins <- unlist(tile(resizePeaks,width=1))

source("~/projects/EPICC/ATAC/4-excluded_samples.R")
excluded_samples <- str_match(excluded_samples, "EPICC_\\s*(.*?)\\s*_C1$")[,2]	

nf_cutsites <- list.files(path = "~/projects/EPICC/ATAC/cutsites",
                          pattern = ".*[1-9]_nucleosome_free-GRCh38-proper_pairs-shifted.bed.gz$",
			                    full.names = TRUE,
			                    recursive=TRUE)

samples <- str_match(nf_cutsites, "EPICC_\\s*(.*?)\\s*_C1_nuc")[,2]
exclude <- samples %in% excluded_samples
samples <- samples[!exclude]
nf_cutsites <- nf_cutsites[!exclude]

combos <- rownames(coldata_region)

# Read in the cutsites
invisible(lapply(1:length(nf_cutsites), function(c){
  gr <- unlist(GRangesList(
    lapply(which(grepl(combos[c],nf_cutsites)), function(s) {
      gr <-import(nf_cutsites[s],
                  format="BED",
                  genome="hg38",
                  which=resizePeaks) 
      gr <- resize(gr, width = 100, fix = "center")
      return(gr)
    })))
  saveRDS(gr,paste0("nf_cov_new/",combos[c],"_peakcov_nf.rds"))
}))

# Get the coverage
cov <- function(libsize,normfacs,dir,name){
  librarysize <- coldata_region$libsize * normfacs
  millions <- librarysize / 1000000
  covs <- lapply(1:length(combos), function(c) {
    cov <- coverage(readRDS(paste0(dir,"_cov_new/", combos[c],"_peakcov_",dir,".rds")))/millions[c]
    cov <- cov[seqlevels(cov) %in% seqlevels(bins)]
    cov <- GenomicRanges::binnedAverage(bins, cov, "score")
    return(cov)})
  names(covs) <- combos
  saveRDS(covs,name)
}
cov("peak_libsize_nf",1,"nf","cpm_peaknorm_cov_nf_1bp.rds")

# get the coverage for region E (normals)
nf_cutsites <- "region_E_nucleosome_free-GRCh38-proper_pairs-shifted_overlap.bed"
librarysize <- 524633802
gr <- import(nf_cutsites,
             format="BED",
             genome="hg38",
             which=resizePeaks) 
gr <- resize(gr, width = 100, fix = "center")
saveRDS(gr,"nf_cov_new/panpatientE_peakcov_nf.rds")
millions <- librarysize / 1000000
cov <- coverage(gr)/millions
cov <- cov[seqlevels(cov) %in% seqlevels(bins)]
cov <- GenomicRanges::binnedAverage(bins, cov, "score")
saveRDS(cov,"panpatientE_cpm_peaknorm_cov_nf_1bp.rds")








