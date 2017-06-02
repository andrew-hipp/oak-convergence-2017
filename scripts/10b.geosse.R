## ML fit, GEOSSE

model.full.4c <- make.geosse(dat.biogeo$tr,
                             dat.biogeo$geosse[tr.spp.4c.discreteClock.noOG$tip.label],
							 sampling.f = c(0.5, 68/154, (13 + 79) / (25 + 89)))
							 # sampled 68 out of estimated 154 Mexican spp, and for North America and the Roburoids, 79/89 and 13/25 respectively
model.4c.no.sAB <- constrain(model.full.4c, sAB ~ 0)
model.4c.noDiff <- constrain(model.full.4c, sA ~ sB, xA ~ xB)

model.4c.p <- starting.point.geosse(dat.biogeo$tr)
model.4c.fitted.full <- find.mle(model.full.4c, model.4c.p)
model.4c.fitted.no.sAB <- find.mle(model.4c.no.sAB, model.4c.p)
model.4c.fitted.noDiff <- find.mle(model.4c.noDiff, model.4c.p)

geosse.4c.results <- round(rbind(full = c(coef(model.4c.fitted.full), lnL = logLik(model.4c.fitted.full), aic = AIC(model.4c.fitted.full)),
            no.sAB = c(coef(model.4c.fitted.no.sAB, T), lnL = logLik(model.4c.fitted.no.sAB), aic = AIC(model.4c.fitted.no.sAB)),
            noDiff = c(coef(model.4c.fitted.noDiff, T), lnL = logLik(model.4c.fitted.noDiff), aic = AIC(model.4c.fitted.noDiff))), 6)
write.csv(geosse.4c.results, 'out/geosse.4c.MLfit.csv')

geosse.4c.anova = anova(model.4c.fitted.full, no.sAB = model.4c.fitted.no.sAB, noDiff = model.4c.fitted.noDiff)

## MCMC, GEOSSE, NO sAB
p.4c.mcmc.no.sAB <- coef(model.4c.fitted.no.sAB)
prior <- make.prior.exponential(1/2)

set.seed(1)
tmp <- mcmc(model.4c.no.sAB, p.4c.mcmc.no.sAB, nsteps=100, prior=prior, w=1, print.every=0)
w <- diff(sapply(tmp[2:7], quantile, c(0.025, 0.975)))
rm(tmp)

mcmc.4c.no.sAB <- mcmc(model.4c.no.sAB, p.4c.mcmc.no.sAB, nsteps = 25000, prior = prior, w = w)

mcmc.4c.diff <- with(mcmc.4c.no.sAB, data.frame(s.diff=sA-sB,
  x.diff=xA-xB, d.diff=dA-dB, div.A=sA-xA, div.B=sB-xB))
colMeans(mcmc.4c.diff > 0)

pdf('out/mcmc.4c.plot.no.sAB.longRun.2017-05-25.pdf')
col1 <- c("red", "orange", "blue", "purple", "black", "gray")
col2 <- col1[c(1,3,5)]
mcmc.4c.diff <- with(mcmc.4c.no.sAB, data.frame(s.diff=sA-sB, x.diff=xA-xB, d.diff=dA-dB))
par(mfrow=c(2,1), mar=c(3, 4, 0, 1))
profiles.plot(mcmc.4c.no.sAB[2:7], col.line=col1, n.br = 100, xlab="", ylab="")
legend("topright", argnames(model.4c.no.sAB), col=col1, lty=1)
profiles.plot(mcmc.4c.diff, col.line=col2, n.br = 100, xlab="", ylab="")
legend("topright", colnames(mcmc.4c.diff), col=col2, lty=1)
title(xlab="rate", ylab="posterior probability density", outer=T, line=-1)
dev.off()

save(mcmc.4c.no.sAB, file = 'out/mcmc.4c.no.sAB.longRun.2017-02-22.Rdata')
