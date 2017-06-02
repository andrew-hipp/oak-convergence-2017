## drought severity mapping

do.rjmcmc = F
writePlots.mcmc = T

if(do.rjmcmc) {
  file.out <- format(Sys.time(), "rjmcmc.%Y-%m-%d")
  tr.im <- drop.tip(tr.spp.4c.discreteClock.noOG,
                    setdiff(tr.spp.4c.discreteClock.noOG$tip.label, row.names(im.means)))
  im.means.tr <- im.means[which(row.names(im.means) %in% tr.im$tip.label), ]
  rjmcmc.bm(tr.im, im.means.tr[, 'Im'], im.means.tr[, 'sem'], ngen=2500000, samp=1000,
            filebase=paste(file.out, 'empirical.sem', sep = '.'),
            simple.start=TRUE, type="jump-rbm")

  rjmcmc.bm(tr.im, im.means.tr[, 'Im'], ngen=1000000, samp=1000,
            filebase=paste(file.out, 'mcmc.sem', sep = '.'),
            simple.start=TRUE, type="jump-rbm")
}

## looking at rjmcmc run

if(writePlots.mcmc) {
  im.rj.emp <- load.rjmcmc("jump-relaxedBM.rjmcmc.2017-05-31.empirical.sem")
  coda::autocorr.plot(im.rj.emp$log, ask=dev.interactive())
  plot(im.rj.emp$log, ask=dev.interactive())

  ## looking at data
  pdf('out/rjmcmc.empiricalErrors.pdf', 12, 8)
  plot(x=im.rj.emp, par="jumps", burnin=0.25, legend=F, show.tip=F, edge.width=3)
  dev.off()
  pdf('out/rjmcmc.empiricalErrors.withLegend.pdf', 12, 8)
  plot(x=im.rj.emp, par="shifts", burnin=0.25, legend=T, show.tip=F, edge.width=3)
  dev.off()


  im.rj.est <- load.rjmcmc("jump-relaxedBM.rjmcmc.2017-05-31.mcmc.sem")
  coda::autocorr.plot(im.rj.est$log, ask=dev.interactive())
  plot(im.rj.est$log, ask=dev.interactive())

  plot(x=im.rj.est, par="jumps", burnin=0.25, legend=TRUE, show.tip=FALSE, edge.width=3)

  im.shifts <- OUshifts(im.means.tr[, 'Im'], tr.im, nmax = 10)
  im.log.shifts <- OUshifts(log(im.means.tr[, 'Im'] - min(im.means.tr[, 'Im']) + 1), tr.im, nmax = 10)
}
