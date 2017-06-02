## various statistics from diversification and niche analyses


## increases in diversification rates
(div.ages.4c.table.pretty['sect. Lobatae, Mexico lineage', 'Diversification rate, epsilon = 0.0'] -
  div.ages.4c.table.pretty['sect. Lobatae', 'Diversification rate, epsilon = 0.0']) /
  div.ages.4c.table.pretty['sect. Lobatae', 'Diversification rate, epsilon = 0.0']

(div.ages.4c.table.pretty['sect. Quercus, Mexico lineage', 'Diversification rate, epsilon = 0.0'] -
  div.ages.4c.table.pretty['sect. Quercus s.s.', 'Diversification rate, epsilon = 0.0']) /
  div.ages.4c.table.pretty['sect. Quercus s.s.', 'Diversification rate, epsilon = 0.0']

  (div.ages.4c.table.pretty['sect. Lobatae, Mexico lineage', 'Diversification rate, epsilon = 0.9'] -
    div.ages.4c.table.pretty['sect. Lobatae', 'Diversification rate, epsilon = 0.9']) /
    div.ages.4c.table.pretty['sect. Lobatae', 'Diversification rate, epsilon = 0.9']

  (div.ages.4c.table.pretty['sect. Quercus, Mexico lineage', 'Diversification rate, epsilon = 0.9'] -
    div.ages.4c.table.pretty['sect. Quercus s.s.', 'Diversification rate, epsilon = 0.9']) /
    div.ages.4c.table.pretty['sect. Quercus s.s.', 'Diversification rate, epsilon = 0.9']

geosse.4c.results

## ordination
eco.mds.fitted <- list(k2 = envfit(eco.mds.multi[[2]]$all, eco.means$all),
                       k5 = envfit(eco.mds.multi[[5]]$all, eco.means$all))


pdf('out/SUPPLEMENT.eco.mds.ordisurf.v2.pdf', 12, 18)
todo <- c('bio10', 'bio11', 'bio4', 'bio12', 'latitude', 'longitude')
layout(matrix(1:6, 3, 2, byrow = TRUE))
eco.mds.surface <- vector('list', length(todo))
for(i in todo) {
  plot(eco.mds.multi[[2]]$all, main = i)
  eco.mds.surface[[i]] <- ordisurf(eco.mds.multi[[2]]$all, eco.means$all[, i], add = TRUE)
  text(-10, 10, paste('r2 = ', round(eco.mds.fitted$k2$vectors$r[i], 3)))
}
dev.off()


pdf('out/eco.mds.fitted5.pdf')
plot(eco.mds.multi[[5]]$all)
plot(eco.mds.fitted$k5)
dev.off()

pdf('out/eco.mds.fitted2.pdf')
plot(eco.mds.multi[[2]]$all)
plot(eco.mds.fitted$k2)
dev.off()

eco.mds.multi[[2]]$all
eco.mds.multi[[5]]$all

## l1ou
eco.shifts.l1ou$convergent
eco.shifts.l1ou.5d.convergent
log(2)  / eco.shifts.l1ou$convergent$alpha # mds1 adaptive half-time
log(2)  / eco.shifts.l1ou.5d.convergent$alpha # mds1 adaptive half-time

## rj.mcmc
pdf('out/rj.mcmc.lnl.pdf')
plot(im.rj.emp$log[, 'lnL'], type = 'l')
dev.off()

mx.wo.branches.rj <- im.rj.emp$phy$hash[which(sapply(paintSubTree(im.rj.emp$phy, fastMRCA(im.rj.emp$phy, 'Quercus_mohriana', 'Quercus_turbinella'), 0)$maps, names) == '0')]
mx.ro.branches.rj <- im.rj.emp$phy$hash[which(sapply(paintSubTree(im.rj.emp$phy, fastMRCA(im.rj.emp$phy, 'Quercus_castanea', 'Quercus_canbyi'), 0)$maps, names) == '0')]
i = c(paste('BASE RATE: ', round(mean(apply(im.rj.emp$rates[-c(1:100), ], 1, mean)), 2) , " (", paste(round(quantile(apply(im.rj.emp$rates[-c(1:100), ], 1, mean), c(0.025, 0.975)), 2), collapse = ",") , ")", sep = ''),
      paste('MX / TX RATE: ', round(mean(apply(im.rj.emp$rates[-c(1:100), c(mx.wo.branches.rj, mx.ro.branches.rj)], 1, mean)), 2) , " (", paste(round(quantile(apply(im.rj.emp$rates[-c(1:100), c(mx.wo.branches.rj, mx.ro.branches.rj)], 1, mean), c(0.025, 0.975)), 2), collapse = ",") , ")", sep = ''),
      paste('MX / TX RATE, white oaks only: ', round(mean(apply(im.rj.emp$rates[-c(1:100), c(mx.wo.branches.rj)], 1, mean)), 2) , " (", paste(round(quantile(apply(im.rj.emp$rates[-c(1:100), c(mx.wo.branches.rj)], 1, mean), c(0.025, 0.975)), 2), collapse = ",") , ")", sep = ''),
      paste('MX / TX RATE, red oaks only: ', round(mean(apply(im.rj.emp$rates[-c(1:100), c(mx.ro.branches.rj)], 1, mean)), 2) , " (", paste(round(quantile(apply(im.rj.emp$rates[-c(1:100), c(mx.ro.branches.rj)], 1, mean), c(0.025, 0.975)), 2), collapse = ",") , ")", sep = ''),
      paste('NON-MX / TX RATE: ', round(mean(apply(im.rj.emp$rates[-c(1:100), setdiff(colnames(im.rj.emp$rates), c(mx.wo.branches.rj, mx.ro.branches.rj))], 1, mean)), 2) , " (", paste(round(quantile(apply(im.rj.emp$rates[-c(1:100), setdiff(colnames(im.rj.emp$rates), c(mx.wo.branches.rj, mx.ro.branches.rj))], 1, mean), c(0.025, 0.975)), 2), collapse = ",") , ")", sep = '')
)
print(i)

## leaf evolution
save(lf.IG10.e.bd.mat,
  lf.gee.e.bd.mat,
  lf.IG10.eb.d.mat,
  lf.gee.eb.d.mat,
  lf.IG10.ed.mat,
  lf.gee.ed.mat,
  lf.fullModels.scaledData,
  file = 'out/lf.logistic.regression.matrices.Rdata')

round(lf.IG10.e.bd.mat,3)
round(lf.gee.e.bd.mat,3)
round(lf.IG10.eb.d.mat,3)
round(lf.gee.eb.d.mat,3)
round(lf.IG10.ed.mat,3)
round(lf.gee.ed.mat,3)

lf.regs.aicw <- cbind(IG10.e.bd = lf.IG10.e.bd.mat[, 'aic.w'],
  gee.e.bd=lf.gee.e.bd.mat[, 'aic.w'],
  IG10.eb.d=lf.IG10.eb.d.mat[, 'aic.w'],
  gee.eb.d=lf.gee.eb.d.mat[, 'aic.w'],
  IG10.ed=lf.IG10.ed.mat[, 'aic.w'],
  gee.ed=lf.gee.ed.mat[, 'aic.w'])

message('\nIves and Garland mean AICw')
apply(lf.regs.aicw[, grep('IG10', dimnames(lf.regs.aicw)[[2]])], 1, mean)
message('\nGEE mean AICw')
apply(lf.regs.aicw[, grep('gee', dimnames(lf.regs.aicw)[[2]])], 1, mean)

lapply(lf.fullModels.scaledData, summary)

sapply(mcmc.2state, function(w) apply(changes(w)$perUnitBrL, 1, function(x) paste(round(x[2], 3), ' (', round(x[1], 3), ",", round(x[3], 3), ')', sep = '')))

ev.out <- c((paste('Evergreen root, transitions to evergreenness: ', mean(mcmc.2.summary$evRoot$count[, '0,1']), " (", paste(quantile(mcmc.2.summary$evRoot$count[, '0,1'], c(0.025,  0.975)), collapse = ','), ")", sep = '')),
             paste('Deciduous root, transitions to evergreenness: ', mean(mcmc.2.summary$decRoot$count[, '0,1']), " (", paste(quantile(mcmc.2.summary$decRoot$count[, '0,1'], c(0.025,  0.975)), collapse = ','), ")", sep = '')
             )
print(ev.out)

wo.mrca <- fastMRCA(tr.leaf, 'Quercus_lobata', 'Quercus_glabrescens')
ro.mrca <- fastMRCA(tr.leaf, 'Quercus_kelloggii', 'Quercus_canbyi')
mx.wo.mrca <- fastMRCA(tr.leaf, 'Quercus_turbinella', 'Quercus_glabrescens')
mx.ro.mrca <- fastMRCA(tr.leaf, 'Quercus_castanea', 'Quercus_canbyi')
mcmc.2.summary$evRoot$ace[as.character(c(wo.mrca, ro.mrca, mx.wo.mrca, mx.ro.mrca)), ]
mcmc.2.summary$decRoot$ace[as.character(c(wo.mrca, ro.mrca, mx.wo.mrca, mx.ro.mrca)), ]
mcmc.2.summary$estimated$ace[as.character(c(wo.mrca, ro.mrca, mx.wo.mrca, mx.ro.mrca)), ]
row.names(lf.ml.sym$lik.anc) <- row.names(lf.ml.bin.sym$lik.anc) <- row.names(mcmc.2.summary$estimated$ace)
lf.ml.sym$lik.anc[as.character(c(wo.mrca, ro.mrca, mx.wo.mrca, mx.ro.mrca)), ]
lf.ml.bin.sym$lik.anc[as.character(c(wo.mrca, ro.mrca, mx.wo.mrca, mx.ro.mrca)), ]
