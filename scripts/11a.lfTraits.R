## updated 2017-03-10
## 2017-05-30: with completed data matrix, the ML reconstruction still makes sense,
##    but MCMC does not. Still need to try binary coding.

ncores <- 14 # set to 1 if you have a windows machine, or lower number of cores depending on available threads

do.ML <- F
do.mcmc <- F
do.summaries <- T

if(do.ML) {
  lf.spp <- intersect(tr.spp.4c.discreteClock.noOG$tip.label, row.names(lf.traits)[which(lf.traits$lfPhenology != '')])
  tr.leaf <- drop.tip(tr.spp.4c.discreteClock.noOG, which(!tr.spp.4c.discreteClock.noOG$tip.label %in% lf.spp))

  lf.ml.sym <- ace(lf1[lf.spp], tr.leaf,
               type = 'discrete',
        			 model = "SYM"
        		   )
  lf.ml.bin.sym <- ace(lf2[lf.spp], tr.leaf,
                  type = 'discrete',
                  model = "SYM"
                  )

  pdf('out/lf.reconstruction.v2-1.panel.pdf', 11, 8.5)
  layout(matrix(1:2, 1))
  for(i in list(lf.ml.sym, lf.ml.bin.sym)) {
    if(dim(i$index.matrix)[1] == 3) pie.colors = structure(c('white', 'grey', 'black'), names = c('Deciduous', 'Brevideciduous', 'Evergreen'))
    if(dim(i$index.matrix)[1] == 2) pie.colors = structure(c('white', 'black'), names = c('Deciduous', 'Evergreen'))
    plot(tr.leaf, show.tip.label = F, edge.width = 0.5)
    nodelabels(pie = i$lik.anc, piecol = pie.colors, cex = 0.8)
    if(dim(i$index.matrix)[1] == 3) tiplabels(pch = 22, bg = pie.colors[lf1[lf.spp] + 1])
    if(dim(i$index.matrix)[1] == 2) tiplabels(pch = 22, bg = pie.colors[lf2[lf.spp] + 1])
    legend(0, 30, names(pie.colors), pch = 22, pt.bg = pie.colors, bty = 'n')
  }
  dev.off()
}

if(do.mcmc) {
  pi.3state = list(evRoot = c(0,0,1),
                   decRoot = c(1,0,0),
                   estimated = 'estimated')
  names(pi.3state$evRoot) <- names(pi.3state$decRoot) <- c(0,1,2)

  pi.2state = list(evRoot = c(0,1),
                   decRoot = c(1,0),
                   estimated = 'estimated')
  names(pi.2state$evRoot) <- names(pi.2state$decRoot) <- c(0,1)

  mcmc.3state <- mclapply(pi.3state, function(thisPi) {
    make.simmap(tr.leaf, x = lf1[lf.spp], nsim=1000, model = "SYM", pi = thisPi)
    },
    mc.cores = ncores)
  names(mcmc.3state) <- names(pi.3state)

  mcmc.2state <- mclapply(pi.2state, function(thisPi) {
    make.simmap(tr.leaf, x = lf2[lf.spp], nsim=1000, model = "SYM", pi = thisPi)
    },
    mc.cores = ncores)
  names(mcmc.2state) <- names(pi.2state)

  mcmc.3.summary <- mclapply(mcmc.3state, describe.simmap, mc.cores = ncores)
  mcmc.2.summary <- mclapply(mcmc.2state, describe.simmap, mc.cores = ncores)

  pdf('out/lf.reconstruction.mcmc.v3-sym.pdf', 12, 16)

  layout(matrix(1:6, 2, 3, byrow = TRUE))

  lapply(mcmc.3.summary, function(x) {
    plot(tr.leaf, show.tip.label = F, edge.width = 0.5)
    nodelabels(pie = x$ace, piecol = c('white', 'grey', 'black'), cex = 0.5)
    tiplabels(pch = 22, bg = c('white', 'grey', 'black')[lf1[lf.spp] + 1], cex = 0.85)
    legend(0, 25, c('Evergreen', 'Brevideciduous', 'Deciduous'), pch = 22, pt.bg = c('black', 'grey', 'white'), bty = 'n')
    })

  lapply(mcmc.2.summary, function(x) {
    plot(tr.leaf, show.tip.label = F, edge.width = 0.5)
    nodelabels(pie = x$ace, piecol = c('white', 'black'), cex = 0.5)
    tiplabels(pch = 22, bg = c('white', 'black')[lf2[lf.spp] + 1], cex = 0.85)
    legend(0, 25, c('Evergreen', 'Deciduous'), pch = 22, pt.bg = c('black', 'white'), bty = 'n')
    })

  dev.off()
}

if(do.summaries) {
  mx.wo.branches <- which(sapply(paintSubTree(tr.leaf, fastMRCA(tr.leaf, 'Quercus_turbinella', 'Quercus_glabrescens'), 0)$maps, names) == '0')
  mx.ro.branches <- which(sapply(paintSubTree(tr.leaf, fastMRCA(tr.leaf, 'Quercus_castanea', 'Quercus_canbyi'), 0)$maps, names) == '0')

  ## redo below to give distributions from MCMC sample
  changes <- function(mcmc, tr = tr.leaf, br.list = list(mx.ro = mx.ro.branches, mx.wo = mx.wo.branches, mx = c(mx.ro.branches, mx.wo.branches),
                                                         not.mx = setdiff(seq(dim(tr.leaf$edge)[1]), c(mx.ro.branches, mx.wo.branches))),
                      quants = c(0.025, 0.50, 0.975)) {
    changeMat <- t(sapply(mcmc, function(x) sapply(x$maps, length) - 1))
    changesRaw <- t(sapply(names(br.list), function(i) quantile(rowSums(changeMat[, br.list[[i]]]), probs = quants)))
    changesPerBranchL <- t(sapply(names(br.list), function(i) quantile(rowSums(changeMat[, br.list[[i]]]) / sum(tr$edge.length[br.list[[i]]]), probs = quants)))
     out = list(raw = changesRaw, perUnitBrL = changesPerBranchL)
     out
   }
}
