
source("horizon_mass_plot_function_commonE_CL6.R")

peaks_promoup <- peaks_original[peaks_original$hypothesis=="gain" & peaks_original$type=="promoter", ]
peaks_promodn <- peaks_original[peaks_original$hypothesis=="loss" & peaks_original$type=="promoter", ]
peaks_enhup <- peaks_original[peaks_original$hypothesis=="gain" & peaks_original$type=="enhancer", ]
peaks_enhdn <- peaks_original[peaks_original$hypothesis=="loss" & peaks_original$type=="enhancer", ]

pdf("~/promo_gain_suppl.pdf", height=5*length(peaks_promoup), width=7*length(patients))
grid.newpage()
pushViewport(viewport(layout = grid.layout(length(peaks_promoup), length(patients))))
for(pat in 1:length(patients)){
for(b in 1:length(peaks_promoup)){
  pushViewport(viewport(layout.pos.col = pat,
                        layout.pos.row = b))
  print(plot_horizon_mass(cov=cov,
                     patient=patients[pat],
                     peak=peaks_promoup[b,]$peak,
                     smooth=40,
                     lwd=1.5,
                     legend=T,
                     type="l"), vp= vplayout(x=pat,y=b))
  popViewport()
}}
dev.off()

pdf("~/promo_loss_suppl.pdf", height=5*length(peaks_promodn), width=7*length(patients))
grid.newpage()
pushViewport(viewport(layout = grid.layout(length(peaks_promodn), length(patients))))
for(pat in 1:length(patients)){
  for(b in 1:length(peaks_promodn)){
    pushViewport(viewport(layout.pos.col = pat,
                          layout.pos.row = b))
    print(plot_horizon(cov=cov,
                       patient=patients[pat],
                       peak=peaks_promodn[b,]$peak,
                       smooth=40,
                       lwd=1.5,
                       legend=T,
                       type="l"), vp= vplayout(x=pat,y=b))
    popViewport()
  }}
dev.off()


pdf("enhancer_gain_suppl.pdf", height=5*length(peaks_enhup), width=7*length(patients))
grid.newpage()
pushViewport(viewport(layout = grid.layout(length(peaks_enhup), length(patients))))
for(pat in 1:length(patients)){
  for(b in 1:length(peaks_enhup)){
    pushViewport(viewport(layout.pos.col = pat,
                          layout.pos.row = b))
    print(plot_horizon_mass(cov=cov,
                       patient=patients[pat],
                       peak=peaks_enhup[b,]$peak,
                       smooth=40,
                       lwd=1.5,
                       legend=T,
                       type="l"), vp= vplayout(x=pat,y=b))
    popViewport()
  }}
dev.off()

pdf("enhancer_loss_suppl.pdf", height=5*length(peaks_enhdn), width=7*length(patients))
grid.newpage()
pushViewport(viewport(layout = grid.layout(length(peaks_enhdn), length(patients))))
for(pat in 1:length(patients)){
  for(b in 1:length(peaks_enhdn)){
    pushViewport(viewport(layout.pos.col = pat,
                          layout.pos.row = b))
    print(plot_horizon_mass(cov=cov,
                       patient=patients[pat],
                       peak=peaks_enhdn[b,]$peak,
                       smooth=40,
                       lwd=1.5,
                       legend=T,
                       type="l"), vp= vplayout(x=pat,y=b))
    popViewport()
  }}
dev.off()

