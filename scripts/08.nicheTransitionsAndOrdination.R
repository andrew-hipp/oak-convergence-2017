## Ecological analyses... rests on the diversification analyses run in '02.treeAnalyses' first to load data

today <- format(Sys.time(), '%Y-%m-%d')
do.mds = TRUE
do.eco.shifts = T

pchSect <- c(21, 25, 22, 22, 23, 24)
colorSect <- c('red', 'white', 'green', 'blue', 'black', 'yellow')
names(colorSect) <- names(pchSect) <- c('Lobatae', 'Protobalanus', 'Quercus - American', 'Quercus - Eurasian', 'Sadlerianae', 'Virentes')

colorRegion <- c('blue', 'black', 'red')
names(colorRegion) <- c("ENA", "MX.SW", "CA")


###################################################

if(do.mds) {
eco.mds <- list(soils = metaMDS(scale(eco.means$soil), 'euclidean', k=2),
                clim = metaMDS(scale(eco.means$clim), 'euclidean', k=2),
				all = metaMDS(scale(cbind(eco.means$soil, eco.means$clim)), 'euclidean', k=2)
        )

eco.mds.multi <- lapply(1:10, function(x) {list(soils = metaMDS(scale(eco.means$soil), 'euclidean', k=x),
                clim = metaMDS(scale(eco.means$clim), 'euclidean', k=x),
        				all = metaMDS(scale(cbind(eco.means$soil, eco.means$clim)), 'euclidean', k=x)
                )
              }
              )
} # close do.mds

eco.mds.stressVec <- sapply(eco.mds.multi, function(x) {sapply(x, function(y) y$stress)})
pdf('out/SUPPL.mds.stress.vec.pdf')
matplot(t(eco.mds.stressVec), type = 'l', col = c('black', 'red', 'blue'))
legend(4, 0.4, row.names(eco.mds.stressVec), col = c('black', 'red', 'blue'))

for(mdsType in c('soils', 'clim', 'all')) {
#  for(mdsDims in list(1:2, c(1,3))) {
    mdsDims = 1:2
    pdf(paste('out/phyloecospace.by.clade', mdsType, paste(mdsDims, collapse = "_"), today, 'pdf', sep = '.'), 8.5, 11)
    layout(matrix(c(1,2), 2, 1))
    wo.points <- which(row.names(eco.mds[[mdsType]]$points) %in% taxa$wo)
    ro.points <- which(row.names(eco.mds[[mdsType]]$points) %in% taxa$ro)
    phylomorphospace(drop.tip(tr.eco, which(!tr.eco$tip.label %in% taxa$wo)), eco.mds[[mdsType]]$points[wo.points, mdsDims], label = 'off', node.size = 0,
                    xlim = range(eco.mds[[mdsType]]$points[c(wo.points, ro.points), mdsDims[1]]),
                    ylim = range(c(-10, eco.mds[[mdsType]]$points[c(wo.points, ro.points), mdsDims[2]])))
    points(eco.mds[[mdsType]]$points[wo.points,mdsDims], pch = 21, col = 'black',
          bg = colorRegion[sect.species.translate[row.names(eco.mds$all$points)[wo.points], 'subclade']], cex = 1.5)
    title(main = 'North American sect Quercus mapped in ecological space')
    legend(4,-5, unique(sect.species.translate[row.names(eco.mds$all$points)[wo.points], 'subclade']), pch = 21, col = 'black',
          pt.bg = colorRegion[unique(sect.species.translate[row.names(eco.mds$all$points)[wo.points], 'subclade'])], bty = 'n')
    phylomorphospace(drop.tip(tr.eco, which(!tr.eco$tip.label %in% taxa$ro)), eco.mds[[mdsType]]$points[ro.points,mdsDims], label = 'off', node.size = 0,
                    xlim = range(eco.mds[[mdsType]]$points[c(wo.points, ro.points), mdsDims[1]]),
                    ylim = range(c(-10, eco.mds[[mdsType]]$points[c(wo.points, ro.points), mdsDims[2]])))
    points(eco.mds[[mdsType]]$points[ro.points,mdsDims], pch = 21, col = 'black',
           bg = colorRegion[sect.species.translate[row.names(eco.mds$all$points)[ro.points], 'subclade']], cex = 1.5)
    title(main = 'North American sect Lobatae mapped in ecological space')
    legend(4,-5, unique(sect.species.translate[row.names(eco.mds$all$points)[ro.points], 'subclade']), pch = 21, col = 'black',
           pt.bg = colorRegion[unique(sect.species.translate[row.names(eco.mds$all$points)[ro.points], 'subclade'])], bty = 'n')

    dev.off()
#}
}

if(do.eco.shifts) {
  eco.shifts <-list (all.1 = OUshifts(eco.mds$all$points[, 1], tr.eco, nmax = 10),
                    all.2 = OUshifts(eco.mds$all$points[, 2], tr.eco, nmax = 10))

  pdf('out/eco.shifts.all1.pdf')
  plot(eco.shifts$all.1, show.tip.label = F)
  dev.off()

  pdf('eco.shifts.all2.pdf')
  plot(eco.shifts$all.2, show.tip.label= F)
  dev.off()


  require(l1ou)
  eco.dat.l1ou <- adjust_data(tr.eco, eco.mds$all$points, normalize = F)
  eco.shifts.l1ou <- list(all.multiVar = estimate_shift_configuration(eco.dat.l1ou$tree, eco.dat.l1ou$Y),
                          all.1 = estimate_shift_configuration(eco.dat.l1ou$tree, eco.dat.l1ou$Y[, 1]),
                          all.2 = estimate_shift_configuration(eco.dat.l1ou$tree, eco.dat.l1ou$Y[, 2]))
  eco.shifts.l1ou$convergent <- estimate_convergent_regimes(eco.shifts.l1ou$all.multiVar)
  #eco.shifts.l1ou$convergent.boots <- l1ou_bootstrap_support(eco.shifts.l1ou$convergent, 100, TRUE, 16, F)

  pdf('out/eco.shifts.l1ou.multiVar.pdf')
  plot(eco.shifts.l1ou$all.multiVar)
  dev.off()

  pdf('out/eco.shifts.l1ou.multiVarConvergent.pdf')
  plot(eco.shifts.l1ou$convergent)
  dev.off()

  eco.dat.l1ou.5d <- adjust_data(tr.eco, eco.mds.multi[[5]]$all$points, normalize = F)
  eco.shifts.l1ou.5d <- estimate_shift_configuration(eco.dat.l1ou.5d$tree, eco.dat.l1ou.5d$Y)
  eco.shifts.l1ou.5d.convergent.2 <- estimate_convergent_regimes(eco.shifts.l1ou.5d, fixed.alpha = F, nCores = 14)
  # fixed.alpha = T used to avoide "system is computationally singular" error
  # https://github.com/khabbazian/l1ou/issues/30

  pdf('out/eco.shifts.l1ou.multiVar.5d.pdf')
  plot(eco.shifts.l1ou.5d)
  dev.off()

  pdf('out/eco.shifts.l1ou.multiVar.5d.Convergent.pdf')
  plot(eco.shifts.l1ou.5d.convergent)
  dev.off()

}
