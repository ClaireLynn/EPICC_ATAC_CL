library(DESeq2)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
options(ggrepel.max.overlaps = Inf,force=1)

patients <- c("C525","C538","C539","C539","C539","C542","C549","C516","C518","C518","C518","C524","C531","C551","C559","C562")

# make comparison names, titles for plot
comparisons <- c("C525"="C525: all region C vs. rest",
		 "C538"="C538: all region D vs. rest",
		 "C539AvCD"="C539: all region A+B1_G3 vs. all region C+D",
		 "C539Avrest"="C539: all region A+B1_G3 vs. rest excluding all region C+D",
		 "C539CDvrest"="C539: all region C+D vs. rest excluding region A+B1_G3",
		 "C542"="C542: all region ABD vs. all region C",
		 "C549"="C549: all region A vs. rest",
		 "C516"="C516: all region B vs. all region A",
		 "C518DvAB"="C518: all region D vs all region A+B",
		 "C518DvC"="C518: all region D vs. all region C",
		 "C518ABvC"="C518: all region A+B vs all region C",
		 "C524"="C524: all region B vs. rest",
		 "C531"="C531: all region B vs. rest",
		 "C551"="C551: rest vs. C1.G4,B1.G3,A1.G6,B1.G7,B1.G2,A1.G9",
		 "C559"="C559: all region B vs. rest",
		 "C562"="C562: res vs. A1_G7")

# make levels for comparisons
comparison_levels <- list("C525"=c("C525.ABD","C525.C"),
		                      "C538"=c("C538.BC","C538.D"),
		                      "C539AvCD"=c("C539.CD","C539.AB1_G3"),
		                      "C539Avrest"=c("C539.B","C539.AB1_G3"),
		                      "C539CDvrest"=c("C539.B","C539.CD"),
		                      "C542"=c("C542.C","C542.ABD"),
                          "C549"=c("C549.BCD","C549.A"),
                          "C516"=c("C516.A","C516.B"),
		                      "C518DvAB"=c("C518.AB","C518.D"),
		                      "C518DvC"=c("C518.C","C518.D"),
		                      "C518ABvC"=c("C518.C","C518.AB"),
                          "C524"=c("C524.ACD","C524.B"),
                          "C531"=c("C531.ACD","C531.B"),
		                      "C551"=c("C551.select","C551.rest"),
                          "C559"=c("C559.CD","C559.B"),
                          "C562"=c("C562.A1_G7","C562.ABCD")
		                       )

# volcano plot function
plot_volcano <- function(res,name,comparison,p,force){        
  volc <- ggplot(data=as.data.frame(res),
  aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(col=Significance),alpha=0.4,size=4) +
  labs(title=comparison, x="Log Fold Change", y ="-log10 Adjusted P.Value")  
 
   volc <- volc + theme_classic(base_size = 24) +
          scale_color_manual(values=c("Down"="cornflowerblue", "NS"="grey", "Up"="brown3"),drop = FALSE) 
  
  volc
  ggsave(filename=paste0(name,"_volc_plot.pdf"),
  plot = last_plot(), scale = 1, width = 25, height = 30, units ="cm")
  volc <- volc + geom_text_repel(inherit.aes = FALSE,force=force,
          data=res[res$padj<p & !is.na(res$new_id),],
  aes(x=log2FoldChange, y=-log10(padj),label = paste0(new_id,"-",type)),
           size=6, segment.color ="grey", segment.alpha=0.5)
  volc
  ggsave(filename=paste0(name,"_volc_plot_labelled.pdf"),
  plot = last_plot(), scale = 1, width = 25, height = 30, units ="cm")
}

count <- "atacseq_counts_per_peak_annot.rds"
model <- as.formula("~purity+sample_type+group")
reduced <- as.formula("~purity+sample_type")
dir <- "effect_of_group_nfreads_LRT"
coldata <- readRDS("coldata_inf.rds")

for (i in 1:length(patients)) {
args <- c(countsf,modelf,reducedf,dirf,patients[i])

# Get coldata
patient <- patients[i]
comparison <- comparisons[i]
comparison_level <- comparison_levels[[i]]
coldat <- coldata[coldata$patient==patient & coldata$group %in% comparison_level,]
# Change levels of the group to ensure direcion of comparison
coldat$group <- factor(coldat$group, levels=comparison_level)

# Get counts
counts <- readRDS(count)
colnames(counts) <- gsub("EPICC_|_C1$","",colnames(counts))
counts <- counts[,rownames(coldat)]

#Only use matching peaks for patient
peak_matches_df<-readRDS("peak_matches_df.rds")
match <- setNames(peak_matches_df[,patient],rownames(peak_matches_df))
match <- names(match[which(match==TRUE)])
counts <- counts[match,]

# Do DESeq
dds <- DESeqDataSetFromMatrix(countData=counts,
         colData=coldat,
         design=model)
dds <- DESeq(dds, test="LRT", reduced=reduced)
res <- results(dds)

pdf(file=paste0(dir,"/",patient,"_",gsub(" ","",as.character(model)[2]),"overlapPeaks_MAplot.pdf"), height=6, width=6)
plotMA(res,alpha=0.05, ylim=c(-2,2.5))
dev.off()

final_reccurence_count <- readRDS("final_reccurence_data.rds")$summary

res <- as.data.frame(res[complete.cases(res),])
res <- merge(final_reccurence_count,res,by.x=1,by.y=0)
res <- res[order(res$padj),]

res$Significance <- "NS"
res$Significance[res$log2FoldChange > 0 & res$padj < 0.05] <- "Up"
res$Significance[res$log2FoldChange < 0 & res$padj < 0.05] <- "Down"

# Write results
saveRDS(res,paste0(dir,"/",names(comparison_levels)[i],"_",gsub(" ","",as.character(model)[2]),"overlapPeaks_deseq_LRT.rds"))

# Plot volcano
plot_volcano(res,name=paste0(dir,"/",names(comparison_levels)[i]),comparison=comparison, p=0.05,force=10)


function(res,name,comparison,p,force){

}
