lf.traits$lf.e.bd <- lf.traits$lf.eb.d <- lf.traits$lf1
lf.traits$lf.e.bd[lf.traits$lf.e.bd == 1] <- 0
lf.traits$lf.eb.d[lf.traits$lf.eb.d == 1] <- 2

lf.traits_lf.ed <- structure(lf.traits$lf1, names = row.names(lf.traits))
lf.traits_lf.ed <- lf.traits_lf.ed[-which(lf.traits_lf.ed == 1)]
lf.traits_lf.ed[lf.traits_lf.ed == 2] <- 1

# lf.traits$lf.ed <- lf.traits$lf.ed[which(lf.traits$lf1 != 1)]
lf.traits$lf.eb.d[lf.traits$lf.eb.d == 2] <- 1
lf.traits$lf.e.bd[lf.traits$lf.e.bd == 2] <- 1
lf.traits$lfPhenology.bin <- lf.traits$lfPhenology
lf.traits$lfPhenology.bin[lf.traits$lfPhenology.bin %in% c('Brevideciduous', "Deciduous, Brevideciduous")] <- 'Deciduous'
#row.names(lf.traits)[row.names(lf.traits) == 'Quercus_michauxi'] <- row.names(im.means)[row.names(im.means) == 'Quercus_michauxi'] <- tr.spp.4c.discreteClock.noOG$tip.label[tr.spp.4c.discreteClock.noOG$tip.label == 'Quercus_michauxi'] <- 'Quercus_michauxii'
#row.names(lf.traits)[row.names(lf.traits) == 'Quercus_magnolifolia'] <- row.names(im.means)[row.names(im.means) == 'Quercus_magnolifolia'] <- tr.spp.4c.discreteClock.noOG$tip.label[tr.spp.4c.discreteClock.noOG$tip.label == 'Quercus_magnfolifolia'] <- 'Quercus_magnoliifolia'
lf.spp <- intersect(tr.spp.4c.discreteClock.noOG$tip.label, row.names(lf.traits)[which(lf.traits$lfPhenology != '')])
# lf.spp.ed <- intersect(tr.spp.4c.discreteClock.noOG$tip.label, row.names(lf.traits)[which(lf.traits$lfPhenology != '')])
lf.spp <- intersect(lf.spp, row.names(eco.means$all))
lf.bin.dat <- cbind(im.means[lf.spp, ], lf.traits[lf.spp, ], eco.means$all[lf.spp, ])
lf.bin.tr <- drop.tip(tr.spp.4c.discreteClock.noOG, which(!tr.spp.4c.discreteClock.noOG$tip.label %in% lf.spp))

lf.bin.dat.rescaled <- cbind(scale(im.means[lf.spp, ]), lf.traits[lf.spp, ], scale(eco.means$all[lf.spp, ]))

lf.spp_ed <- intersect(lf.spp, names(lf.traits_lf.ed))
lf.bin.dat_ed <- as.data.frame(cbind(im.means[lf.spp_ed,], eco.means$all[lf.spp_ed, ], lf.ed = lf.traits_lf.ed[lf.spp_ed]))
lf.bin.tr_ed <- drop.tip(tr.spp.4c.discreteClock.noOG, which(!tr.spp.4c.discreteClock.noOG$tip.label %in% lf.spp_ed))

lf.bin.e.bd.gee <- list(
  e.bd.11.10 = phyloglm(lf.e.bd ~ bio11 + bio10, lf.bin.dat, phy=lf.bin.tr, btol = 35),
  e.bd.im.11.10 = phyloglm(lf.e.bd ~ Im + bio11 + bio10, lf.bin.dat, phy=lf.bin.tr, btol = 35),
  e.bd.im.11.4.15 = phyloglm(lf.e.bd ~ Im + bio11 + bio4 + bio15, lf.bin.dat, phy=lf.bin.tr, btol = 35),
  e.bd.im.11.4.10.15 = phyloglm(lf.e.bd ~ Im + bio11 + bio4 + bio10 + bio15, lf.bin.dat, phy=lf.bin.tr, btol = 35),
  e.bd.im.10.interaction = phyloglm(lf.e.bd ~ Im + bio10 + Im*bio10, lf.bin.dat, lf.bin.tr, btol = 35),
  e.bd.im.11.interaction = phyloglm(lf.e.bd ~ Im + bio11 + Im*bio11, lf.bin.dat, lf.bin.tr, btol = 35)
#  e.bd.im.11.10.scaled = phyloglm(lf.e.bd ~ Im + bio11 + bio10, lf.bin.dat.rescaled, phy=lf.bin.tr, btol = 35),
#  e.bd.im.11.4.15.scaled = phyloglm(lf.e.bd ~ Im + bio11 + bio4 + bio15, lf.bin.dat.rescaled, phy=lf.bin.tr, btol = 35),
#  e.bd.im.11.4.10.15.scaled = phyloglm(lf.e.bd ~ Im + bio11 + bio4 + bio10 + bio15, lf.bin.dat.rescaled, phy=lf.bin.tr, btol = 35),
#  e.bd.im.10.interaction.scaled = phyloglm(lf.e.bd ~ Im + bio10 + Im*bio10, lf.bin.dat.rescaled, lf.bin.tr, btol = 35)
)

lf.bin.e.bd.IG10 <- list(
  e.bd.11.10 = phyloglm(lf.e.bd ~ bio11 + bio10, lf.bin.dat, phy=lf.bin.tr, btol = 35, method = 'logistic_IG10'),
  e.bd.im.11.10 = phyloglm(lf.e.bd ~ Im + bio11 + bio10, lf.bin.dat, phy=lf.bin.tr, btol = 35, method = 'logistic_IG10'),
  e.bd.im.11.4.15 = phyloglm(lf.e.bd ~ Im + bio11 + bio4 + bio15, lf.bin.dat, phy=lf.bin.tr, btol = 35, method = 'logistic_IG10'),
  e.bd.im.11.4.10.15 = phyloglm(lf.e.bd ~ Im + bio11 + bio4 + bio10 + bio15, lf.bin.dat, phy=lf.bin.tr, btol = 35, method = 'logistic_IG10'),
  e.bd.im.10.interaction = phyloglm(lf.e.bd ~ Im + bio10 + Im*bio10, lf.bin.dat, lf.bin.tr, btol = 35, method = 'logistic_IG10'),
  e.bd.im.11.interaction = phyloglm(lf.e.bd ~ Im + bio11 + Im*bio11, lf.bin.dat, lf.bin.tr, btol = 35, method = 'logistic_IG10')
#  e.bd.im.11.10.scaled = phyloglm(lf.e.bd ~ Im + bio11 + bio10, lf.bin.dat.rescaled, phy=lf.bin.tr, btol = 35, method = 'logistic_IG10'),
#  e.bd.im.11.4.15.scaled = phyloglm(lf.e.bd ~ Im + bio11 + bio4 + bio15, lf.bin.dat.rescaled, phy=lf.bin.tr, btol = 35, method = 'logistic_IG10'),
#  e.bd.im.11.4.10.15.scaled = phyloglm(lf.e.bd ~ Im + bio11 + bio4 + bio10 + bio15, lf.bin.dat.rescaled, phy=lf.bin.tr, btol = 35, method = 'logistic_IG10'),
#  e.bd.im.10.interaction.scaled = phyloglm(lf.e.bd ~ Im + bio10 + Im*bio10, lf.bin.dat.rescaled, lf.bin.tr, btol = 35, method = 'logistic_IG10')
)

lf.bin.eb.d.gee <- list(
  eb.d.11.10 = phyloglm(lf.eb.d ~  bio11 + bio10, lf.bin.dat, phy=lf.bin.tr, btol = 35),
  eb.d.im.11.10 = phyloglm(lf.eb.d ~ Im + bio11 + bio10, lf.bin.dat, phy=lf.bin.tr, btol = 35),
  eb.d.im.11.4.15 = phyloglm(lf.eb.d ~ Im + bio11 + bio4 + bio15, lf.bin.dat, phy=lf.bin.tr, btol = 35),
  eb.d.im.11.4.10.15 = phyloglm(lf.eb.d ~ Im + bio11 + bio4 + bio10 + bio15, lf.bin.dat, phy=lf.bin.tr, btol = 35),
  eb.d.im.10.interaction = phyloglm(lf.eb.d ~ Im + bio10 + Im*bio10, lf.bin.dat, lf.bin.tr, btol = 35),
  eb.d.im.11.interaction = phyloglm(lf.eb.d ~ Im + bio11 + Im*bio11, lf.bin.dat, lf.bin.tr, btol = 35)
#  eb.d.im.11.10.scaled = phyloglm(lf.eb.d ~ Im + bio11 + bio10, lf.bin.dat.rescaled, phy=lf.bin.tr, btol = 35),
#  eb.d.im.11.4.15.scaled = phyloglm(lf.eb.d ~ Im + bio11 + bio4 + bio15, lf.bin.dat.rescaled, phy=lf.bin.tr, btol = 35),
#  eb.d.im.11.4.10.15.scaled = phyloglm(lf.eb.d ~ Im + bio11 + bio4 + bio10 + bio15, lf.bin.dat.rescaled, phy=lf.bin.tr, btol = 35),
#  eb.d.im.10.interaction.scaled = phyloglm(lf.eb.d ~ Im + bio10 + Im*bio10, lf.bin.dat.rescaled, lf.bin.tr, btol = 35)
)

lf.bin.eb.d.IG10 <- list(
  eb.d.11.10 = phyloglm(lf.eb.d ~ bio11 + bio10, lf.bin.dat, phy=lf.bin.tr, btol = 35, method = 'logistic_IG10'),
  eb.d.im.11.10 = phyloglm(lf.eb.d ~ Im + bio11 + bio10, lf.bin.dat, phy=lf.bin.tr, btol = 35, method = 'logistic_IG10'),
  eb.d.im.11.4.15 = phyloglm(lf.eb.d ~ Im + bio11 + bio4 + bio15, lf.bin.dat, phy=lf.bin.tr, btol = 35, method = 'logistic_IG10'),
  eb.d.im.11.4.10.15 = phyloglm(lf.eb.d ~ Im + bio11 + bio4 + bio10 + bio15, lf.bin.dat, phy=lf.bin.tr, btol = 35, method = 'logistic_IG10'),
  eb.d.im.10.interaction = phyloglm(lf.eb.d ~ Im + bio10 + Im*bio10, lf.bin.dat, lf.bin.tr, btol = 35, method = 'logistic_IG10'),
  eb.d.im.11.interaction = phyloglm(lf.eb.d ~ Im + bio11 + Im*bio11, lf.bin.dat, lf.bin.tr, btol = 35, method = 'logistic_IG10')
#  eb.d.im.11.10.scaled = phyloglm(lf.eb.d ~ Im + bio11 + bio10, lf.bin.dat.rescaled, phy=lf.bin.tr, btol = 35, method = 'logistic_IG10'),
#  eb.d.im.11.4.15.scaled = phyloglm(lf.eb.d ~ Im + bio11 + bio4 + bio15, lf.bin.dat.rescaled, phy=lf.bin.tr, btol = 35, method = 'logistic_IG10'),
#  eb.d.im.11.4.10.15.scaled = phyloglm(lf.eb.d ~ Im + bio11 + bio4 + bio10 + bio15, lf.bin.dat.rescaled, phy=lf.bin.tr, btol = 45, method = 'logistic_IG10'),
#  eb.d.im.10.interaction.scaled = phyloglm(lf.eb.d ~ Im + bio10 + Im*bio10, lf.bin.dat.rescaled, lf.bin.tr, btol = 35, method = 'logistic_IG10')
)

lf.bin.ed.gee <- list(
  ed.11.10 = phyloglm(lf.ed ~ bio11 + bio10, lf.bin.dat_ed, phy=lf.bin.tr_ed, btol = 35),
  ed.im.11.10 = phyloglm(lf.ed ~ Im + bio11 + bio10, lf.bin.dat_ed, phy=lf.bin.tr_ed, btol = 35),
  ed.im.11.4.15 = phyloglm(lf.ed ~ Im + bio11 + bio4 + bio15, lf.bin.dat_ed, phy=lf.bin.tr_ed, btol = 35),
  ed.im.11.4.10.15 = phyloglm(lf.ed ~ Im + bio11 + bio4 + bio10 + bio15, lf.bin.dat_ed, phy=lf.bin.tr_ed, btol = 35),
  ed.im.10.interaction = phyloglm(lf.ed ~ Im + bio10 + Im*bio10, lf.bin.dat_ed, lf.bin.tr_ed, btol = 35),
  ed.im.11.interaction = phyloglm(lf.ed ~ Im + bio11 + Im*bio11, lf.bin.dat_ed, lf.bin.tr_ed, btol = 35)
#  ed.im.11.10.scaled = phyloglm(lf.ed ~ Im + bio11 + bio10, lf.bin.dat_ed.rescaled, phy=lf.bin.tr_ed, btol = 45),
#  ed.im.11.4.15.scaled = phyloglm(lf.ed ~ Im + bio11 + bio4 + bio15, lf.bin.dat_ed.rescaled, phy=lf.bin.tr_ed, btol = 45),
#  ed.im.11.4.10.15.scaled = phyloglm(lf.ed ~ Im + bio11 + bio4 + bio10 + bio15, lf.bin.dat_ed.rescaled, phy=lf.bin.tr_ed, btol = 45),
#  ed.im.10.interaction.scaled = phyloglm(lf.ed ~ Im + bio10 + Im*bio10, lf.bin.dat_ed.rescaled, lf.bin.tr_ed, btol = 45)
)

lf.bin.ed.IG10 <- list(
  ed.11.10 = phyloglm(lf.ed ~ bio11 + bio10, lf.bin.dat_ed, phy=lf.bin.tr_ed, btol = 35, method = 'logistic_IG10'),
  ed.im.11.10 = phyloglm(lf.ed ~ Im + bio11 + bio10, lf.bin.dat_ed, phy=lf.bin.tr_ed, btol = 35, method = 'logistic_IG10'),
  ed.im.11.4.15 = phyloglm(lf.ed ~ Im + bio11 + bio4 + bio15, lf.bin.dat_ed, phy=lf.bin.tr_ed, btol = 35, method = 'logistic_IG10'),
  ed.im.11.4.10.15 = phyloglm(lf.ed ~ Im + bio11 + bio4 + bio10 + bio15, lf.bin.dat_ed, phy=lf.bin.tr_ed, btol = 35, method = 'logistic_IG10'),
  ed.im.10.interaction = phyloglm(lf.ed ~ Im + bio10 + Im*bio10, lf.bin.dat_ed, lf.bin.tr_ed, btol = 35, method = 'logistic_IG10'),
  ed.im.11.interaction = phyloglm(lf.ed ~ Im + bio11 + Im*bio11, lf.bin.dat_ed, lf.bin.tr_ed, btol = 35, method = 'logistic_IG10')
#  ed.im.11.10.scaled = phyloglm(lf.ed ~ Im + bio11 + bio10, lf.bin.dat_ed.rescaled, phy=lf.bin.tr_ed, btol = 30, method = 'logistic_IG10'),
#  ed.im.11.4.15.scaled = phyloglm(lf.ed ~ Im + bio11 + bio4 + bio15, lf.bin.dat_ed.rescaled, phy=lf.bin.tr_ed, btol = 30, method = 'logistic_IG10'),
#  ed.im.11.4.10.15.scaled = phyloglm(lf.ed ~ Im + bio11 + bio4 + bio10 + bio15, lf.bin.dat_ed.rescaled, phy=lf.bin.tr_ed, btol = 45),
#  ed.im.10.interaction.scaled = phyloglm(lf.ed ~ Im + bio10 + Im*bio10, lf.bin.dat_ed.rescaled, lf.bin.tr_ed, btol = 30, method = 'logistic_IG10')
)

lf.bin.e.bd.mple <- list(
  e.bd.11.10 = phyloglm(lf.e.bd ~ bio11 + bio10, lf.bin.dat, phy=lf.bin.tr, btol = 35, method = 'logistic_MPLE'),
  e.bd.im.11.10 = phyloglm(lf.e.bd ~ Im + bio11 + bio10, lf.bin.dat, phy=lf.bin.tr, btol = 35, method = 'logistic_MPLE'),
  e.bd.im.11.4.15 = phyloglm(lf.e.bd ~ Im + bio11 + bio4 + bio15, lf.bin.dat, phy=lf.bin.tr, btol = 35, method = 'logistic_MPLE'),
  e.bd.im.11.4.10.15 = phyloglm(lf.e.bd ~ Im + bio11 + bio4 + bio10 + bio15, lf.bin.dat, phy=lf.bin.tr, btol = 35, method = 'logistic_MPLE'),
  e.bd.im.10.interaction = phyloglm(lf.e.bd ~ Im + bio10 + Im*bio10, lf.bin.dat, lf.bin.tr, btol = 35, method = 'logistic_MPLE'),
  e.bd.im.11.interaction = phyloglm(lf.e.bd ~ Im + bio11 + Im*bio11, lf.bin.dat, lf.bin.tr, btol = 35, method = 'logistic_MPLE')
)

lf.IG10.e.bd.mat <- lf.gee.e.bd.mat <- lf.mple.e.bd.mat <- matrix(NA, nrow = length(lf.bin.e.bd.IG10), ncol = 9,
                    dimnames = list(names(lf.bin.e.bd.IG10), c('Im', 'bio4', 'bio11', 'bio10', 'bio15', 'Im:bio10', 'Im:bio11', 'aic', 'aic.w')))
lf.IG10.eb.d.mat <- lf.gee.eb.d.mat <- matrix(NA, nrow = length(lf.bin.eb.d.IG10), ncol = 9,
                    dimnames = list(names(lf.bin.eb.d.IG10), c('Im', 'bio4', 'bio11', 'bio10', 'bio15', 'Im:bio10', 'Im:bio11', 'aic', 'aic.w')))
lf.IG10.ed.mat <- lf.gee.ed.mat <- matrix(NA, nrow = length(lf.bin.ed.IG10), ncol = 9,
                    dimnames = list(names(lf.bin.ed.IG10), c('Im', 'bio4', 'bio11', 'bio10', 'bio15', 'Im:bio10', 'Im:bio11', 'aic', 'aic.w')))

for(i in names(lf.bin.e.bd.IG10)) {
  p.ig = summary(lf.bin.e.bd.IG10[[i]])$coefficients[-c(1), 'p.value']
  lf.IG10.e.bd.mat[i, names(p.ig)] <- p.ig
  p.gee = summary(lf.bin.e.bd.gee[[i]])$coefficients[-c(1), 'p.value']
  lf.gee.e.bd.mat[i, names(p.gee)] <- p.gee
  p.mple = summary(lf.bin.e.bd.mple[[i]])$coefficients[-c(1), 'p.value']
  lf.mple.e.bd.mat[i, names(p.gee)] <- p.gee

}
for(i in names(lf.bin.eb.d.IG10)) {
  message(i)
  p.ig = summary(lf.bin.eb.d.IG10[[i]])$coefficients[-c(1), 'p.value']
  lf.IG10.eb.d.mat[i, names(p.ig)] <- p.ig
  p.gee = summary(lf.bin.eb.d.gee[[i]])$coefficients[-c(1), 'p.value']
  lf.gee.eb.d.mat[i, names(p.gee)] <- p.gee
}
for(i in names(lf.bin.ed.IG10)) {
  p.ig = summary(lf.bin.ed.IG10[[i]])$coefficients[-c(1), 'p.value']
  lf.IG10.ed.mat[i, names(p.ig)] <- p.ig
  p.gee = summary(lf.bin.ed.gee[[i]])$coefficients[-c(1), 'p.value']
  lf.gee.ed.mat[i, names(p.gee)] <- p.gee
}

lf.IG10.e.bd.mat[, 'aic'] <- sapply(lf.bin.e.bd.IG10, AIC)
lf.gee.e.bd.mat[, 'aic'] <- sapply(lf.bin.e.bd.gee, AIC)
lf.mple.e.bd.mat[, 'aic'] <- sapply(lf.bin.e.bd.gee, AIC)
lf.IG10.eb.d.mat[, 'aic'] <- sapply(lf.bin.eb.d.IG10, AIC)
lf.gee.eb.d.mat[, 'aic'] <- sapply(lf.bin.eb.d.gee, AIC)
lf.IG10.ed.mat[, 'aic'] <- sapply(lf.bin.ed.IG10, AIC)
lf.gee.ed.mat[, 'aic'] <- sapply(lf.bin.ed.gee, AIC)


lf.IG10.e.bd.mat[, 'aic.w'] <- exp(-0.5 * lf.IG10.e.bd.mat[, 'aic'] - min (lf.IG10.e.bd.mat[, 'aic'])) / sum(exp(-0.5 * lf.IG10.e.bd.mat[, 'aic'] - min (lf.IG10.e.bd.mat[, 'aic'])))
lf.gee.e.bd.mat[, 'aic.w'] <- exp(-0.5 * lf.gee.e.bd.mat[, 'aic'] - min (lf.gee.e.bd.mat[, 'aic'])) / sum(exp(-0.5 * lf.gee.e.bd.mat[, 'aic'] - min (lf.gee.e.bd.mat[, 'aic'])))
lf.mple.e.bd.mat[, 'aic.w'] <- exp(-0.5 * lf.mple.e.bd.mat[, 'aic'] - min (lf.mple.e.bd.mat[, 'aic'])) / sum(exp(-0.5 * lf.mple.e.bd.mat[, 'aic'] - min (lf.mple.e.bd.mat[, 'aic'])))
lf.IG10.eb.d.mat[, 'aic.w'] <- exp(-0.5 * lf.IG10.eb.d.mat[, 'aic'] - min (lf.IG10.eb.d.mat[, 'aic'])) / sum(exp(-0.5 * lf.IG10.eb.d.mat[, 'aic'] - min (lf.IG10.eb.d.mat[, 'aic'])))
lf.gee.eb.d.mat[, 'aic.w'] <- exp(-0.5 * lf.gee.eb.d.mat[, 'aic'] - min (lf.gee.eb.d.mat[, 'aic'])) / sum(exp(-0.5 * lf.gee.eb.d.mat[, 'aic'] - min (lf.gee.eb.d.mat[, 'aic'])))
lf.IG10.ed.mat[, 'aic.w'] <- exp(-0.5 * lf.IG10.ed.mat[, 'aic'] - min (lf.IG10.ed.mat[, 'aic'])) / sum(exp(-0.5 * lf.IG10.ed.mat[, 'aic'] - min (lf.IG10.ed.mat[, 'aic'])))
lf.gee.ed.mat[, 'aic.w'] <- exp(-0.5 * lf.gee.ed.mat[, 'aic'] - min (lf.gee.ed.mat[, 'aic'])) / sum(exp(-0.5 * lf.gee.ed.mat[, 'aic'] - min (lf.gee.ed.mat[, 'aic'])))

lf.fullModels.scaledData <- list(
  e.bd.IG10.im.11.4.10.15 = phyloglm(lf.e.bd ~ Im + bio11 + bio4 + bio10 + bio15, lf.bin.dat.rescaled, phy=lf.bin.tr, btol = 35, method = 'logistic_IG10'),
  eb.d.IG10.im.11.4.10.15 = phyloglm(lf.eb.d ~ Im + bio11 + bio4 + bio10 + bio15, lf.bin.dat.rescaled, phy=lf.bin.tr, btol = 35, method = 'logistic_IG10'),
  e.bd.gee.im.11.4.10.15 = phyloglm(lf.e.bd ~ Im + bio11 + bio4 + bio10 + bio15, lf.bin.dat.rescaled, phy=lf.bin.tr, btol = 35),
  eb.d.gee.im.11.4.10.15 = phyloglm(lf.eb.d ~ Im + bio11 + bio4 + bio10 + bio15, lf.bin.dat.rescaled, phy=lf.bin.tr, btol = 35),
  ed.IG10.im.11.4.10.15 = phyloglm(lf.eb.d ~ Im + bio11 + bio4 + bio10 + bio15, lf.bin.dat.rescaled[-grep('Brevideciduous', lf.bin.dat.rescaled$lfPhenology), ], phy=drop.tip(lf.bin.tr, row.names(lf.bin.dat.rescaled)[grep('Brevideciduous', lf.bin.dat.rescaled$lfPhenology)]), btol = 35, method = 'logistic_IG10'),
  ed.gee.im.11.4.10.15 = phyloglm(lf.eb.d ~ Im + bio11 + bio4 + bio10 + bio15, lf.bin.dat.rescaled[-grep('Brevideciduous', lf.bin.dat.rescaled$lfPhenology), ], phy=drop.tip(lf.bin.tr, row.names(lf.bin.dat.rescaled)[grep('Brevideciduous', lf.bin.dat.rescaled$lfPhenology)]), btol = 35)
  )

write.csv(lf.IG10.e.bd.mat, 'out/lf.IG10.e.bd.mat.csv')

pdf('out/leafPhenology.boxplots.pdf', 15, 8)
layout(matrix(1:5, 1, 5))
boxplot(Im ~ lfPhenology, lf.bin.dat, main = 'Im ~ lfPhen')
boxplot(bio11 ~ lfPhenology, lf.bin.dat, main = 'bio11 ~ lfPhen')
boxplot(bio10 ~ lfPhenology, lf.bin.dat, main = 'bio10 ~ lfPhen')
boxplot(bio4 ~ lfPhenology, lf.bin.dat, main = 'bio4 ~ lfPhen')
boxplot(bio15 ~ lfPhenology, lf.bin.dat, main = 'bio15 ~ lfPhen')
dev.off()

pdf('out/leafPhenology.binary.boxplots.pdf', 15, 8)
layout(matrix(1:5, 1, 5))
boxplot(Im ~ lfPhenology.bin, lf.bin.dat, main = 'Im ~ lfPhen')
boxplot(bio11 ~ lfPhenology.bin, lf.bin.dat, main = 'bio11 ~ lfPhen')
boxplot(bio10 ~ lfPhenology.bin, lf.bin.dat, main = 'bio10 ~ lfPhen')
boxplot(bio4 ~ lfPhenology.bin, lf.bin.dat, main = 'bio4 ~ lfPhen')
boxplot(bio15 ~ lfPhenology.bin, lf.bin.dat, main = 'bio15 ~ lfPhen')
dev.off()

lf.bin.out.eb.d <- phylolm(lf.eb.d ~ Im + bio11, lf.bin.dat, lf.bin.tr)
#sort(im.m)
