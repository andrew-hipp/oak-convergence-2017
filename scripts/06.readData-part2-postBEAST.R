## A second round of data import after:
##  - BEAST analysis of jackknife trees,
##  - Fossil calibrations
##  - Generation of consensus tree with error bars in treeannotator
##  - Running cleaned specimen data through Infomap Bioregions

## tree reading and cleaning
tr.spp.4c.discreteClock.beast = phyloch:::read.beast('out/trs.calib.jackknife.4.annotated.tre')
tr.spp.2c.discreteClock.beast = phyloch:::read.beast('out/trs.calib.jackknife.2.annotated.tre')
tr.spp.4c.discreteClock = read.nexus('out/trs.calib.jackknife.4.annotated.tre')
tr.spp.2c.discreteClock = read.nexus('out/trs.calib.jackknife.2.annotated.tre')
tr.spp.4c.discreteClock.noOG <- drop.tip(tr.spp.4c.discreteClock, which(!tr.spp.4c.discreteClock$tip.label %in% taxa$ingroup))
tr.spp.2c.discreteClock.noOG <- drop.tip(tr.spp.2c.discreteClock, which(!tr.spp.2c.discreteClock$tip.label %in% taxa$ingroup))

## bioregions reading and cleaning
bioregions <- read.table('../data/bioregions.out.v2/bioregions.table.2017-02-22_presence-absence.edit.txt', sep = '|', as.is = T, header = T, row.names = 1)
bioregions$X <- NULL
bioregions.coords <- read.table('../data/bioregions.out.v2/bioregions.table.2017-02-22_bioregions-coords.txt', as.is  = T)
bioregions.coords <- bioregions.coords[, c(2,1)]
bioregions.trans <- read.table('../data/bioregions.out.v2/bioregions.trans.tsv', as.is = T, header = T, row.names = 1)
names(bioregions) <- bioregions.trans$region

dat.biogeo <- list(tr = drop.tip(tr.spp.4c.discreteClock.noOG, which(!tr.spp.4c.discreteClock.noOG$tip.label %in% row.names(tip.geog))),
                     dat = tip.geog[which(row.names(tip.geog) %in% tr.spp.4c.discreteClock.noOG$tip.label), ])
dat.biogeo$geosse <- dat.biogeo$dat[, 'M'] + apply(dat.biogeo$dat[, c('E', 'C', 'U')], 1, function(x) ifelse(any(x == 1), 2, 0))
dat.biogeo$geosse[dat.biogeo$geosse == 3] <- 0
names(dat.biogeo$geosse) <- row.names(dat.biogeo$dat)
