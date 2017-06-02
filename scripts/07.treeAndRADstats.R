## tree and dataset stats

readRad = F
makeMaps = TRUE

if(readRad) {
rad <- read.pyRAD('../data/rad/c85d6m20p3.full.loci')
}

tr.main <- tr[[1]]
tr.main$tip.label <- gsub(" ", "_", tr.main$tip.label)
rad.reads <- read.delim('../data/tables/sequence.counts.2017-01-17.tsv',as.is = T, row.names = 1)
row.names(rad.reads) <- paste('>', gsub('.fq.gz|.fq|_R1', '', row.names(rad.reads)), sep = '')
fna.spp <- readLines('../data/tables/fna.oaks.txt')

dat.renamed$spName <- label.elements(dat.renamed$cleanedSpecies.spOnly, delim = " |_", returnNum = 2)
sppToDo <- sort(unique(dat.renamed$spName[grep('Quercus', dat.renamed$cleanedSpecies)]))
N <- length(sppToDo)

message('FNA taxa missing from tree:')
print(setdiff(fna.spp, label.elements(tr.main$tip.label, delim = '[_|]', returnNum = 1:2, returnDelim = ' ') ))

message('Relative support for alternative models:')
print(apply(chronosTest.mat, 2, mean))

## tree stats:
quercus.tips <- grep("Quercus", tr.main$tip.label, value = TRUE)
ingroup.tips <- grep('libani|acutissima|trojana|baronii|Litho|Castan', tr.main$tip.label, value = T, invert = T)
message(paste('total taxa', length(sapply(strsplit(tr.main$tip.label, '|', fixed = T), function(x) x[1])), sep = ' - '))
message(paste('outgroup tips', paste(grep('Quercus', tr.main$tip.label, invert = TRUE, value = TRUE), collapse = ', '), sep = ' - '))
message(paste('tips, ingroup', length(sapply(strsplit(ingroup.tips, '_', fixed = T), function(x) x[2])), sep = ' - '))
message(paste('spp, ingroup', length(unique(sapply(strsplit(ingroup.tips, '_', fixed = T), function(x) x[2]))), sep = ' - '))
message(paste('mean tips, ingroup', mean(sort(table(sapply(strsplit(ingroup.tips,"_", fixed = T), function(x) x[2])))), sep  = ' - ')) ## mean number of tips / ingroup sp
message(paste('median tips, ingroup', mean(sort(table(sapply(strsplit(ingroup.tips,"_", fixed = T), function(x) x[2])))), sep  = ' - ')) ## mean number of tips / ingroup sp
message(paste('sd tips, ingroup', sd(sort(table(sapply(strsplit(ingroup.tips,"_", fixed = T), function(x) x[2])))), sep = ' - ')) ## sd, number of tips / ingroup sp
message(paste('spp with exactly 2 tips, ingroup', sum(sort(table(sapply(strsplit(ingroup.tips,"_", fixed = T), function(x) x[2])))==2), sep = ' - ')) ## number of ingroup spp with two tips
message(paste('spp with exactly 3 tips, ingroup', sum(sort(table(sapply(strsplit(ingroup.tips,"_", fixed = T), function(x) x[2])))==3), sep = ' - ')) ## number of ingroup spp with two tips
message(paste('spp with exactly 4 tips, ingroup', sum(sort(table(sapply(strsplit(ingroup.tips,"_", fixed = T), function(x) x[2])))==4), sep = ' - ')) ## number of ingroup spp with two tips
message(paste('spp with exactly 5 tips, ingroup', sum(sort(table(sapply(strsplit(ingroup.tips,"_", fixed = T), function(x) x[2])))==5), sep = ' - ')) ## number of ingroup spp with two tips
message(paste('spp with > 5 tips, ingroup', sum(sort(table(sapply(strsplit(ingroup.tips,"_", fixed = T), function(x) x[2])))>5), sep = ' - ')) ## number of ingroup spp with two tips
message('Species with highest sample sizes:')
print(tail(sort(table(sapply(strsplit(ingroup.tips,"_", fixed = T), function(x) x[2])))), 10)

tdi <- summary.by.elements(tr.main, fixed = FALSE, delim = '[_|]', returnNum = 2)$disparity.mat
message(paste('monophyletic spp:', dim(tdi[which(tdi[, 'disparity'] == 0 & tdi[, 'count'] > 1), ])[1]) )# rows are number of monophyletic spp
message(paste('nonmonophyletic spp:', dim(tdi[which(tdi[, 'disparity'] > 0 & tdi[, 'count'] > 1), ])[1]) )# rows are non-monophyletic spp

summary.latLong <- data.frame(tips = sapply(sppToDo, function(x) sum(dat.renamed$spName == x)),
                          geoN = integer(N),
                          geoMax = numeric(N),
                          geoMin = numeric(N),
                          geoMean = numeric(N),
                          row.names = sppToDo)

if(makeMaps) {
    dir.create('out/species.maps')
    for(i in sppToDo[which(summary.latLong$tips > 1)]) {
      temp.lat <- as.numeric(dat.renamed[grep(i, dat.renamed$cleanedSpecies), 'latitude_georef'])
      temp.long <- as.numeric(dat.renamed[grep(i, dat.renamed$cleanedSpecies), 'longitude_georef'])
      temp.lat.long <- cbind(Long = temp.long[!is.na(temp.long)],
                             Lat = temp.lat[!is.na(temp.lat)])
      summary.latLong[i, 'geoN'] <- dim(temp.lat.long)[1]
      if(summary.latLong[i, 'geoN'] < 2) next
      else {
        temp.results <- haversine(temp.lat.long)
        summary.latLong[i, c('geoMax', 'geoMin', 'geoMean')] <- c(max(temp.results), min(temp.results), mean(temp.results))
        pdf(paste('out/species.maps/', i, '.pdf', sep = ''), 11, 8.5)
        map('world')
        points(temp.lat.long, pch = 19, cex = 1, col = 'blue')
        text(-120, -50, paste('N =', dim(temp.lat.long)[1], '\nMax D =',  round(summary.latLong[i, 'geoMax'],1), 'km'), cex = 1.5)
        dev.off()
      }
    }
}

## sampling distances:
summary.distances <- numeric(0)
summary.distances.names <- character(0)

tdi.v.km.dat <- cbind(summary.latLong[intersect(row.names(tdi), row.names(summary.latLong)[which(summary.latLong$tips > 1)]),],
                      tdi[intersect(row.names(tdi), row.names(summary.latLong)[which(summary.latLong$tips > 1)]),])

cor.test(tdi.v.km.dat$geoMax, tdi.v.km.dat$disparity, method = 'spearman')
cor.test(tdi.v.km.dat[which(tdi.v.km.dat$disparity < 10), 'geoMax'],
         tdi.v.km.dat[which(tdi.v.km.dat$disparity < 10), 'disparity'],
         method = 'spearman')
cor.test(tdi.v.km.dat[which(row.names(tdi.v.km.dat) != 'glabrescens'), 'geoMax'],
         tdi.v.km.dat[which(row.names(tdi.v.km.dat) != 'glabrescens'), 'disparity'],
         method = 'spearman')


message('\nDISTANCE STATS')
geoN.2 <- summary.latLong$geoN > 1
message(paste('number of sp with more than one georeferenced site', sum(geoN.2), sep = ' - '))
message(paste('average maximum distance between spp with more than one georeferenced site',
              mean(summary.latLong$geoMax[which(geoN.2)]), sep = ' - '))
message(paste('median maximum distance between spp with more than one georeferenced site',
              median(summary.latLong$geoMax[which(geoN.2)]), sep = ' - '))
message(paste('proportion of spp with more than one georeferenced site for which maximum distance > 100 km',
              round(sum(summary.latLong$geoMax[which(geoN.2)] > 100) / sum(geoN.2), 2), sep = ' - '))
message(paste('proportion of spp with more than one georeferenced site for which maximum distance <10 km',
              round(sum(summary.latLong$geoMax[which(geoN.2)] < 10) / sum(geoN.2), 2), sep = ' - '))


## rad stats:
message('\nRAD STATS\n')
range(rad.reads$reads)
mean(rad.reads$reads)
sd(rad.reads$reads)
print(rads.full) # gives summary
range(rads.full$radSummary$locus.lengths)
mean(rads.full$radSummary$locus.lengths)
sd(rads.full$radSummary$locus.lengths)

loc.nums <- sort(rowSums(rads.full$radSummary$inds.mat[row.names(rad.reads)[rad.reads$includedFinalTrees], ]))
range(loc.nums)
min(loc.nums)/dim(rads.full$radSummary$inds.mat)[2]
loc.nums[2]/dim(rads.full$radSummary$inds.mat)[2]
