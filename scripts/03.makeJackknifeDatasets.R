## making jackknife samples for constraint tree
## 2016-06-08
## updated 2017-01-03 to use final trees

trPath <- ("/home/andrew/Dropbox/NSF-OAKS-1146488-core-files/Papers/2017-AmOaksFull/ANALYSES/analyses.submission.2017-05-17/data/phylogeny/rellSearch/")
trName <- ('bipart.oakTree.2017.tre.uniqueSpp.originalTips.tre')
tr.spp.origTips <- read.tree(paste(trPath, trName, sep = ''))

if(!all(c('rads.full.mat', 'rads.ingroup.mat', 'rads.noPonticae.mat') %in% ls())) load('../data/rad/allLociFiles.RData') # loads the clustering results from 2016-01-29: c85d6m20p3

loci.m20 <- which(colSums(rads.full$radSummary$inds.mat[tr.spp.origTips$tip.label, ]) >= 20)
halfLoci <- as.integer(length(loci.m20) / 2)

toDo = 1:100

if(!'jackknife.m20' %in% dir(trPath)) {
  dir.create(paste(trPath, 'jackknife.m20', sep = ''))
  dir.create(paste(trPath, 'jackknife.m20.loci', sep = ''))
}

file.copy(paste(trPath, trName, sep = ''), paste(trPath, 'jackknife.m20/', trName, sep = ''))

for(i in toDo) {
  message(paste('making jackknife dataset', i))
  locSample <- sample(colnames(rads.full$radSummary$inds.mat), halfLoci, replace = FALSE)
  rad2phy(rads.full.mat, inds = tr.spp.origTips$tip.label,
        loci = locSample,
        outfile = paste(trPath, 'jackknife.m20/jackknifeDataset.', i, '.phy', sep = ''))
  message('just finished rad2phy, trying to write loci')
  writeLines(locSample, paste(trPath, 'jackknife.m20.loci/jackknifeDataset.', i, '.locSample.txt', sep = ''))
  message('just finished writing loci, back around to next i')
}

dir.create(paste(trPath, 'jackknife.m20/out', sep = ''))

writeLines(paste('raxmlHPC-PTHREADS-AVX -f e -T 7 -p 12345 -t ',
                 trName, ' ',
                 '-m GTRGAMMA -s jackknifeDataset.', toDo, '.phy ',
                 '-n jackknifeTree.', toDo, '.tre ',
                 '-w ', trPath, 'jackknife.m20/out', sep = ''), paste(trPath,'jackknife.m20/raxmlOptimTest.v2.sh', sep = ''))
write.tree(tr.spp.origTips, paste(trPath, 'jackknife.m20/', trName, sep = ''))
writeLines(paste('rm jackknifeDataset.', toDo, '.phy.reduced ; gzip jackknifeDataset.', toDo, '.phy', sep = ''), paste(trPath, 'jackknife.m20/raxmlFileCleanup.v2.sh', sep = ''))
