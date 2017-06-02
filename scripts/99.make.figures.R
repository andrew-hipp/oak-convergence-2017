## make figures
## MAKE SURE TO RUN THIS AFTER ALL OTHER SCRIPTS HAVE BEEN RUN THROUGH 11c
## for New Phyt, 3.5" wide for one column, 7" wide for full page (it seems), and 9.5" max height

path = format(Sys.time(), 'figures.%Y-%m-%d.%H%M')
dir.create(path)
regionColors.figs = c(CA = 'deepskyblue', ENA = 'maroon', MX.SW = 'orange', U = 'black')

# Fig. 1. Sampling map
source('../scripts/99.mapMaking.v2.R')

# Fig. 2. Maximum likelihood phylogeny of the American oak clade, with major biogeographic transitions, leaf morphology, and divergence time estimates. Numbers above the branches are non-parametric bootstrap support values estimated using the fast bootstrapping approach in RAxML.
source('../scripts/99.drawTree.v3.R')

# Fig. 3. Speciation, extinction, and dispersal rates for lineages of Mexico / Central America and American lineages outside of Mexico and North America, estimated under the Geographic State Speciation and Extinction model (GeoSSE) model.
geosse.colors <- c(regionColors.figs[c('MX.SW', 'ENA')], 'blue', 'purple', 'black', 'gray')
pdf(paste(path, '/FIG03.geosse.v3.pdf', sep = ''), 3.5, 3.5)
par(mar=c(3,3,0,0) + 0.1)
profiles.plot(mcmc.4c.no.sAB[2:7], col.line=geosse.colors,
              n.br = 150,
              #xlab="rate", ylab="Posterior probability density",
              cex.lab = 0.5, cex.axis = 0.5,
              lwd = c(2), lty = c(1),
              opacity = 0.5
              )
legend("topright", c('Speciation rate, Mexico', 'Speciation rate, non-Mexico',
                     'Extinction rate, Mexico', 'Extinction rate, non-Mexico',
                     'Dispersal from Mexico northward', 'Dispersal to Mexico'),
                     col=geosse.colors,
                     lty = 1, lwd = 2, bty = 'n', cex = 0.7)
title(xlab="Rate", ylab="Posterior probability density", line=-1, outer = TRUE, cex.lab = 0.7)
dev.off()

# Fig. 4. Ordination of taxa in ecological space, combined soils and climatic data, with phylogeny overlaid. A. Red oaks. B. White oaks.
pdf(paste(path, '/FIG04.ordination.all.pdf', sep = ''), 3.5, 7)
par(mfrow=c(2,1),
    mar = c(1,0,0,0) + 0.1,
    oma = c(2,3,0,0) + 0.1)
wo.points <- which(row.names(eco.mds[["all"]]$points) %in% taxa$wo)
ro.points <- which(row.names(eco.mds[["all"]]$points) %in% taxa$ro)

ordisurf(eco.mds.multi[[2]]$all, eco.means$all[, 'bio11']/10,
        lwd.cl = 0.5, col = 'gray80', main = '', cex = 0, add = FALSE, axes = F)
phylomorphospace(drop.tip(tr.eco, which(!tr.eco$tip.label %in% taxa$wo)), eco.mds[["all"]]$points[wo.points, mdsDims],
                label = 'off', node.size = 0, lwd = 0.75, xlab = '', ylab = '', axes = F,
                xlim = range(eco.mds[["all"]]$points[c(wo.points, ro.points), mdsDims[1]]),
                ylim = range(c(-10, eco.mds[["all"]]$points[c(wo.points, ro.points), mdsDims[2]])), add = TRUE)
box()
axis(side = 2, cex.axis = 0.5)
points(eco.mds[["all"]]$points[wo.points,mdsDims], pch = 21, col = 'black',
      bg = regionColors.figs[sect.species.translate[row.names(eco.mds$all$points)[wo.points], 'subclade']], cex = 1)
text(-10.5,10, substitute(paste("section ", italic('Quercus'))), cex = 0.6, pos = 4)
#text(-11, 8.75, 'Contours = mean temp, coldest quarter', cex = 0.5, pos = 4)
# title(main = 'North American sect Quercus mapped in ecological space')

ordisurf(eco.mds.multi[[2]]$all, eco.means$all[, 'bio11']/10,
        lwd.cl = 0.5, col = 'gray80', main = '', cex = 0, add = FALSE, axes = F)
box()
phylomorphospace(drop.tip(tr.eco, which(!tr.eco$tip.label %in% taxa$ro)), eco.mds[["all"]]$points[ro.points,mdsDims],
                label = 'off', node.size = 0, lwd = 0.75, xlab = '', ylab = '', axes = F,
                xlim = range(eco.mds[["all"]]$points[c(wo.points, ro.points), mdsDims[1]]),
                ylim = range(c(-10, eco.mds[["all"]]$points[c(wo.points, ro.points), mdsDims[2]])),
                add = TRUE)
axis(side = 1, cex.axis = 0.5)
axis(side = 2, cex.axis = 0.5)
points(eco.mds[["all"]]$points[ro.points,mdsDims], pch = 21, col = 'black',
       bg = regionColors.figs[sect.species.translate[row.names(eco.mds$all$points)[ro.points], 'subclade']], cex = 1)
text(-10.5,10, substitute(paste("section ", italic('Lobatae'))), cex = 0.6, pos = 4)
#text(-11, 8.75, 'Contours = latitude', cex = 0.5, pos = 4)
# title(main = 'North American sect Lobatae mapped in ecological space')
legend("bottomright", c("C", "E", "M"),
       pch = 21, col = 'black',
       pt.bg = regionColors.figs[c('CA', 'ENA', 'MX.SW')],
       bty = 'n', cex = 0.7)
title(ylab="MDS axis 2", outer=T, line=2)
title(xlab="MDS axis 1", outer=T, line=1)
dev.off()

# Fig. 5. Adaptive transitions in ecological space as estimated from the two axes of the K=2 NMDS ordination. Changes in color on the tree signify transitions in ecological space that exceed expectations under a Brownian motion model of character evolution. Analysis using the phylogenetic lasso (implemented in the l1ou) detected 17 transitions in ecological space, and subsequent analysis of these transitions in an adaptive framework identified nine selective regimes, including convergent shifts in the two Mexican clades (white oaks and red oaks).
pdf(paste(path, '/FIG05.eco.shifts.l1ou.multiVarConvergent.pdf', sep = ''), 7, 9)
par(mar = c(0,0,0,0) + 0.1)
colors.conv <- c(sample(rainbow(length(eco.shifts.l1ou$convergent$shift.configuration))), "maroon")
mx.regs = which(names(eco.shifts.l1ou$convergent$shift.configuration) == '1')
colors.conv[mx.regs] <- 'orange'
plot(eco.shifts.l1ou$convergent, edge.shift.ann = FALSE, asterisk = FALSE, palette = colors.conv)
lastP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tr.mrca <- mrca(eco.shifts.l1ou$convergent$tree)
shN <- c("Mexico Lobatae" = tr.mrca['Quercus_canbyi', 'Quercus_emoryi'],
            "Mexico Quercus" = tr.mrca['Quercus_mohriana', 'Quercus_turbinella'],
            "Virentes" = tr.mrca['Quercus_geminata', 'Quercus_fusiformis'],
            "mostly Texas" = tr.mrca['Quercus_mohriana', 'Quercus_vaseyana'],
            "mostly Arizona" = tr.mrca['Quercus_arizonica', 'Quercus_turbinella'],
            "Protobalanus" = tr.mrca['Quercus_palmeri', 'Quercus_tomentella']
            )
for(i in names(shN)) text(lastP$xx[shN[i]], lastP$yy[shN[i]],
                          i, adj = c(1.1,-0.5),
                          offset = c(0,1),
                          cex = 0.7,
                          font = 4)
dev.off()

# Table 1. Divergence time estimates with four fossil constraints, estimated clade diversity, and absolute diversification rates calculated using Magallon and Sanderson’s (2001) estimator.
write.csv(div.ages.4c.table.pretty, paste(path, '/TABLE01.diversification.ages.4c.csv', sep = ''))

# Fig. S1. Infomap Bioregions map of biogeographic regions.
pdf(paste(path, '/FIG.S01.infomaps.pdf', sep = ''), 7, 6)

a <- ggplot()
a = a + geom_polygon(data = regions.df, aes(long, lat, group = group, fill = factor(bioregio)))
a = a + geom_path(data = regions.df, aes(long, lat, group = group), colour = 'gray50', size = 0.2)
a = a + geom_map(data = ourArea, map = ourArea, aes(x = long, y = lat, map_id = region), colour = 'black', fill = NA, size = 0.5)
a = a + geom_map(data = ourArea.states, map = ourArea.states, aes(x = long, y = lat, map_id = region), colour = 'black', fill = NA, size = 0.5)
a = a + xlim(-125, -63) + ylim(8, 49)
a = a + xlab("Longitude") + ylab('Latitude')
a = a + scale_fill_manual("Bioregions",
                           labels = bioregions.trans$region,
                           values = c("#babdd1", "#d6b874", "#f0a6a3", "#93ba9b", "#b3946d", "#a7b27c", "#c79fab", "#9aadab", "#c9a394", "#cc8c79", "#e3b089")
                           )

a <- a + theme(legend.position = c(0.90, 0.30), legend.background = element_rect(fill=NA))
print(a)
dev.off()

# Fig. S2. Maximum likelihood phylogeny of all tips sampled, with untransformed branch lengths.
pdf(paste(path, '/FIGS02.allTips.pdf', sep = ''), 8.5, 11)
par(mar = rep(0, 4))
plot(tr, cex = 0.2, show.node.label = TRUE)
dev.off()

# Fig. S3. Plot of non-metric multidimensional scaling stress against ordination dimensions (K).
pdf(paste(path, '/FIGS03.mds.stress.vec.pdf', sep = ''))
matplot(t(eco.mds.stressVec), type = 'l', lwd = 2, lty = 'solid', col = c('black', 'red', 'blue'),
        xlab = "K (number of ordination axes)", ylab = "NMDS stress")
legend(7, 0.4, c('Soils data', 'BIOCLIM data', 'Combined data'), lwd = 2, lty = 'solid', bty = 'n', col = c('black', 'red', 'blue'))
dev.off()

# Fig. S4. Surface plots, raw ordination data plotted over NMDS plot. This might go to main text, pruned to four variables (lat, long, bio11, bio4) and geographic regions of taxa colored.
pdf(paste(path, '/FIGS04.eco.mds.ordisurf.v2.pdf', sep = ''), 9.5, 9.5)
todo <- c('bio11', 'bio4', 'latitude', 'longitude')
layout(matrix(1:4, 2, 2, byrow = TRUE))
eco.mds.surface <- vector('list', length(todo))
for(i in todo) {
  plot(eco.mds.multi[[2]]$all, main = i)
  eco.mds.surface[[i]] <- ordisurf(eco.mds.multi[[2]]$all, eco.means$all[, i], add = TRUE)
  text(-10, 10, paste('r2 = ', round(eco.mds.fitted$k2$vectors$r[i], 3)), cex = 0.5)
}
dev.off()

# Fig. S5. NMDS plots for soils and climate separately.
pdf(paste(path, '/FIGS05.soil.and.clim.ordinations.pdf', sep = ''), 8.5, 11)
layout(matrix(1:4, 2, 2))
for(mdsType in c('soils', 'clim')) {
#  for(mdsDims in list(1:2, c(1,3))) {
    mdsDims = 1:2
    wo.points <- which(row.names(eco.mds[[mdsType]]$points) %in% taxa$wo)
    ro.points <- which(row.names(eco.mds[[mdsType]]$points) %in% taxa$ro)
    phylomorphospace(drop.tip(tr.eco, which(!tr.eco$tip.label %in% taxa$wo)), eco.mds[[mdsType]]$points[wo.points, mdsDims], label = 'off', node.size = 0,
                    xlim = range(eco.mds[[mdsType]]$points[c(wo.points, ro.points), mdsDims[1]]),
                    ylim = range(c(-10, eco.mds[[mdsType]]$points[c(wo.points, ro.points), mdsDims[2]])))
    points(eco.mds[[mdsType]]$points[wo.points,mdsDims], pch = 21, col = 'black',
          bg = colorRegion[sect.species.translate[row.names(eco.mds$all$points)[wo.points], 'subclade']], cex = 1.5)
    title(main = paste('sect Quercus,', mdsType,  'space'))
    legend(2.5,-7, unique(sect.species.translate[row.names(eco.mds$all$points)[wo.points], 'subclade']), pch = 21, col = 'black',
          pt.bg = colorRegion[unique(sect.species.translate[row.names(eco.mds$all$points)[wo.points], 'subclade'])], bty = 'n')
    phylomorphospace(drop.tip(tr.eco, which(!tr.eco$tip.label %in% taxa$ro)), eco.mds[[mdsType]]$points[ro.points,mdsDims], label = 'off', node.size = 0,
                    xlim = range(eco.mds[[mdsType]]$points[c(wo.points, ro.points), mdsDims[1]]),
                    ylim = range(c(-10, eco.mds[[mdlastP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tr.mrca <- mrca(eco.shifts.l1ou$convergent$tree)
shN <- c("Mexico Lobatae" = tr.mrca['Quercus_canbyi', 'Quercus_emoryi'],
            "Mexico Quercus" = tr.mrca['Quercus_mohriana', 'Quercus_turbinella'],
            "Virentes" = tr.mrca['Quercus_geminata', 'Quercus_fusiformis'],
            "mostly Texas" = tr.mrca['Quercus_mohriana', 'Quercus_vaseyana'],
            "mostly Arizona" = tr.mrca['Quercus_arizonica', 'Quercus_turbinella'],
            "Protobalanus" = tr.mrca['Quercus_palmeri', 'Quercus_tomentella']
            )
for(i in names(shN)) text(lastP$xx[shN[i]], lastP$yy[shN[i]],
                          i, adj = c(1.1,-0.5),
                          offset = c(0,1),
                          cex = 0.7,
                          font = 4)
sType]]$points[c(wo.points, ro.points), mdsDims[2]])))
    points(eco.mds[[mdsType]]$poedge.shift.ann = FALSE, asterisk = FALSE, ints[ro.points,mdsDims], pch = 21, col = 'black',
           bg = colorRegion[sect.species.translate[row.names(eco.mds$all$points)[ro.points], 'subclade']], cex = 1.5)
    title(main = paste('sect Lobatae,', mdsType,  'space'))
    legend(2.5,-7, unique(sect.species.translate[row.names(eco.mds$all$points)[ro.points], 'subclade']), pch = 21, col = 'black',
           pt.bg = colorRegion[unique(sect.species.translate[row.names(eco.mds$all$points)[ro.points], 'subclade'])], bty = 'n')
}
dev.off()
#}

# Fig. S6. Ornstein-Uhlenbeck adaptive transitions on the K = 5 NMDS ordination axes.
pdf(paste(path, '/FIGS06.eco.shifts.l1ou.multiVar.5d.Convergent.pdf', sep = ''), 8.5, 9)
plot(eco.shifts.l1ou.5d.convergent, edge.shift.ann = FALSE, asterisk = FALSE, cex = 0.4)
dev.off()

# Fig. S7. Model-averaged rate of evolution along moisture gradient (Im), reconstructed using reversible-jump Markov chain Monte Carlo (rjMCMC).
## TO LABEL NODES, TREE AND LEGEND WERE MERGED IN ADOBE ILLUSTRATOR.
pdf(paste(path, '/FIGS07.rjmcmc.empErrors.treeOnly.pdf', sep = ''), 8.5, 9)
plot(x=im.rj.emp, par="shifts", burnin=0.25, legend=F, show.tip=T, edge.width=3, cex = 0.5)

lastP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tr.mrca <- mrca(im.rj.emp$phy)
shN <- c("Mexico Lobatae" = tr.mrca['Quercus_canbyi', 'Quercus_emoryi'],
            "Mexico Quercus" = tr.mrca['Quercus_mohriana', 'Quercus_turbinella'],
            "Virentes" = tr.mrca['Quercus_geminata', 'Quercus_fusiformis'],
            "mostly Texas" = tr.mrca['Quercus_mohriana', 'Quercus_vaseyana'],
            "mostly Arizona" = tr.mrca['Quercus_arizonica', 'Quercus_turbinella'],
            "Protobalanus" = tr.mrca['Quercus_palmeri', 'Quercus_tomentella']
            )
for(i in names(shN)) text(lastP$xx[shN[i]], lastP$yy[shN[i]],
                          i, adj = c(1.1,-0.5),
                          offset = c(0,1),
                          cex = 0.7,
                          font = 4)

dev.off()

pdf(paste(path, '/FIGS07.rjmcmc.empErrors.legend.pdf', sep = ''), 8.5, 9)
plot(x=im.rj.emp, par="shifts", burnin=0.25, legend=T, show.tip=T, edge.width=3, cex = 0.5)
dev.off()

# Fig. S8. Stochastic mapping of leaf phenology.
pdf(paste(path, '/FIGS08.lf.reconstruction.mcmc.v4.pdf', sep = ''), 8.5, 9)
par(mar = c(0,0,0,0))
tr.mrca <- mrca(tr.leaf)
layout(matrix(1:4, 1, 4, byrow = TRUE))
lapply(mcmc.2.summary, function(x) {
  plot(tr.leaf, show.tip.label = F, edge.width = 0.5)
  nodelabels(pie = x$ace, piecol = c('white', 'black'), cex = 0.8, lwd = 0.4)
  tiplabels(pch = 22, bg = c('white', 'black')[lf2[tr.leaf$tip.label] + 1], cex = 1, lwd = 0.4)
  legend(0, 8, c('Evergreen', 'Deciduous'), pch = 22, pt.bg = c('black', 'white'), bty = 'n', cex = 0.8)

  lastP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  shN <- c("Mexico Lobatae" = tr.mrca['Quercus_canbyi', 'Quercus_emoryi'],
              "Mexico Quercus" = tr.mrca['Quercus_mohriana', 'Quercus_turbinella'],
              "Virentes" = tr.mrca['Quercus_geminata', 'Quercus_fusiformis'],
              "mostly Texas" = tr.mrca['Quercus_mohriana', 'Quercus_vaseyana'],
              "mostly Arizona" = tr.mrca['Quercus_arizonica', 'Quercus_turbinella'],
              "Protobalanus" = tr.mrca['Quercus_palmeri', 'Quercus_tomentella']
              )
  for(i in names(shN)) text(lastP$xx[shN[i]], lastP$yy[shN[i]],
                            i, adj = c(1.1,-0.5),
                            offset = c(0,1),
                            cex = 0.7,
                            font = 3)

  if(x$ace[1,1] == 0) text(0, tr.leaf$Nnode, '(a)', font = 2)
  if(x$ace[1,1] == 1) text(0, tr.leaf$Nnode, '(b)', font = 2)
  if(!x$ace[1,1] %in% c(0,1)) text(0, tr.leaf$Nnode, '(c)', font = 2)
  })

  pie.colors = structure(c('white', 'grey', 'black'), names = c('Deciduous', 'Brevideciduous', 'Evergreen'))
  plot(tr.leaf, show.tip.label = F, edge.width = 0.5)
  nodelabels(pie = lf.ml.sym$lik.anc, piecol = pie.colors, cex = 0.8, lwd = 0.4)
  tiplabels(pch = 22, bg = pie.colors[lf1[tr.leaf$tip.label] + 1], lwd = 0.4)
  legend(0, 8, rev(names(pie.colors)), pch = 22, pt.bg = rev(pie.colors), bty = 'n', cex = 0.8)
  lastP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  shN <- c("Mexico Lobatae" = tr.mrca['Quercus_canbyi', 'Quercus_emoryi'],
              "Mexico Quercus" = tr.mrca['Quercus_mohriana', 'Quercus_turbinella'],
              "Virentes" = tr.mrca['Quercus_geminata', 'Quercus_fusiformis'],
              "mostly Texas" = tr.mrca['Quercus_mohriana', 'Quercus_vaseyana'],
              "mostly Arizona" = tr.mrca['Quercus_arizonica', 'Quercus_turbinella'],
              "Protobalanus" = tr.mrca['Quercus_palmeri', 'Quercus_tomentella']
              )
  for(i in names(shN)) text(lastP$xx[shN[i]], lastP$yy[shN[i]],
                            i, adj = c(1.1,-0.5),
                            offset = c(0,1),
                            cex = 0.7,
                            font = 3)
  text(0, tr.leaf$Nnode, '(d)', font = 2)

dev.off()

# Fig. S9. Leaf phenology correlates with selected climatic variables.
pdf(paste(path, '/FIGS09.leafPhenology.boxplots.pdf', sep = ''), 8.5, 11)
par(mfrow = c(5,4))
lf.bin.dat.temp <- lf.bin.dat
lf.bin.dat.temp$lfPhenology[lf.bin.dat.temp$lfPhenology == "Deciduous"] <- 'D'
lf.bin.dat.temp$lfPhenology[lf.bin.dat.temp$lfPhenology == "Brevideciduous"] <- 'B'
lf.bin.dat.temp$lfPhenology[lf.bin.dat.temp$lfPhenology == "Evergreen"] <- 'E'
lapply(c('Im', paste('bio', 1:19, sep = '')), function(x)
  boxplot(formula(paste(x, '~ lfPhenology')), lf.bin.dat.temp, main = x)
  )
rm(lf.bin.dat.temp)
dev.off()

# Table S1. Samples included in this study.
tr.tips <- tr.uniques[[1]]$sppOnlyOrigTip
tr.tips <- gsub(">", "", tr.tips)
tr.tips <- gsub(".barcodeStripped", "", tr.tips)
tr.tips <- gsub(".nameFixed", "", tr.tips)
tr.tips <- gsub(".techRep", "", tr.tips)

dat$inSingleTipsTree <- dat$tip.label %in% tr.tips
write.csv(dat, paste(path, '/TableS01.samples.csv', sep = ''))

# Table S2. Specimen records included in this study, including specimens excluded during data cleanup.
write.csv(moisture.index, paste(path, '/TableS02.specimenRecords.csv', sep = ''))

# Table S3. Species list, with sectional classification, biogeographic coding, and leaf phenology coding.
tableS3 <- cbind(section = sect.species.translate[tr.spp.4c.discreteClock$tip.label, 'Section'],
                 tip.geog[tr.spp.4c.discreteClock$tip.label, ],
                 lf.traits[tr.spp.4c.discreteClock$tip.label, c('lfPhenology', 'lf.eb.d', 'lf.e.bd')])
write.csv(tableS3, paste(path, '/TableS03.spp.lf.biogeo.sect.csv', sep = ''))

# Table S4. Divergence time estimates with two fossil constraints (cf. Table 1).
write.csv(div.ages.2c.table.pretty, paste(path, '/TableS04.diversification.ages.2c.csv', sep = ''))

# Table S5. AIC weights for logistic regression models
write.csv(lf.regs.aicw, paste(path, '/TableS05.regressionModelAICc.csv', sep = ''))
# Table S6. Coefficients for logistic regression models
write.csv(coef(summary(lf.fullModels.scaledData$e.bd.IG10.im.11.4.10.15)),
          paste(path, '/TableS06.e.bd.IG10.csv', sep = ''))
