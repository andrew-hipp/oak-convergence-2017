## making final datasets
## this is run in the background before phylogenetic analysis is conducted, and generates
##   the datasets used in this paper. In the end, only the 'rads.full.mat' was reported on,
##   as sensitivity analyses based on the other matrices had no effect. For a disucssion
##   of the effects of taxon removal on topology of the white oaks, see McVay, Hipp and Manos
##   2017, https://doi.org/10.1098/rspb.2017.0300
## 2016-03-01

library(RADami)

if(!all(c('rads.full.mat', 'rads.ingroup.mat', 'rads.noPonticae.mat') %in% ls())) load('../data/rad/allLociFiles.RData') # loads the clustering results from 2016-01-29: c85d6m20p3
removals <- readLines('../data/rad/removals.2016.11.30') # these are individuals to be removed for final analysis

rad2phy(rads.full.mat, inds = row.names(rads.full.mat)[-which(row.names(rads.full.mat) %in% removals)],
        loci = names(which(colSums(rads.full$radSummary$inds.mat) > 19)), outfile = 'allTaxa.c85d6m20p3.2016.12.08.phy')
rad2phy(rads.ingroup.mat, inds = row.names(rads.ingroup.mat)[-which(row.names(rads.ingroup.mat) %in% removals)],
        loci = names(which(colSums(rads.ingroup$radSummary$inds.mat) > 19)), outfile = 'ingroup.OOA.c85d6m20p3.2016.12.08.phy')
rad2phy(rads.noPonticae.mat, inds = row.names(rads.noPonticae.mat)[-which(row.names(rads.noPonticae.mat) %in% removals)],
        loci = names(which(colSums(rads.noPonticae$radSummary$inds.mat) > 19)), outfile = 'noPonticae.OOA.c85d6m20p3.2016.12.08.phy')
