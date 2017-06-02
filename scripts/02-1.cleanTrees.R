## updated 2016-06-14 to add in a tree that is unique to spp, but has the original tips on it
## removed interactives because this is now part of paper pipeline
## updated 2017-01-03 to :
##   1. use just a single tree and update to our final taxon list
##   2. no longer add a date stamp to the files it writes

trPath <- ('../data/phylogeny/rellSearch/')
trName <- c('bipart.oakTree.2017.tre') ## always rename the tree you want to this name, and set the path above

ogRenamer <- c('>DM46', '>DM71', '>DM74') ## these are the outgroups
# dat <- read.delim('../data/tables/plate.names.2016.11.30.tsv', as.is = T)

## read data
tr <- read.tree(paste(trPath, trName[1], sep = ''))
if(class(tr) == "phylo") {
  tr <- list(tr)
  class(tr) <- 'multiPhylo'
  }

## do some renaming and export trees
for(trNum in 1:length(tr)) {
  tr[[trNum]] <- ladderize(root(tr[[trNum]], ogRenamer))

	tips.orig <- tr[[trNum]]$tip.label

  tr[[trNum]]$tip.label <- gsub(">", "", tr[[trNum]]$tip.label)
	tr[[trNum]]$tip.label <- gsub(".barcodeStripped", "", tr[[trNum]]$tip.label)
	tr[[trNum]]$tip.label <- gsub(".nameFixed", "", tr[[trNum]]$tip.label)
	tr[[trNum]]$tip.label <- gsub(".techRep", "", tr[[trNum]]$tip.label)

  tips.clean <- tr$tip.label

	tipsIndex <- match(tr[[trNum]]$tip.label, dat$tip.label)
  tip.trans <- structure(paste(dat$cleanedSpecies[tipsIndex], tr[[trNum]]$tip.label, dat$geoPol[tipsIndex], sep = ' | '), names = tr[[trNum]]$tip.label)
	dat.renamed <- dat[tipsIndex, ]
	row.names(dat.renamed) <- as.character(tip.trans)
	tr[[trNum]]$tip.label <- as.character(tip.trans)
  dat.renamed$tips.orig <- tips.orig
	}

write.tree(tr, paste(trPath, trName, '.allTips.tre', sep = ''))

makeUniques <- function(tr, dat) {
  tr.uniques <- list(all = drop.tip(tr, tr$tip.label[duplicated(dat$cleanedSpecies[tipsIndex])]))
	tr.uniques$taxaOnly <- tr.uniques$all
	tr.uniques$taxaOnly$tip.label <- sapply(tr.uniques$all$tip.label, function(x) strsplit(x, " | ", fixed = T)[[1]][1])
	#tr.uniques$sppOnlyOrigTips <- tr.uniques$sppOnly <- tr.uniques$taxaOnly
	# tr.uniques$sppOnly <- tr.uniques$all
  sppOnly.vector <- sapply(tr.uniques$taxaOnly$tip.label, function(x) paste(strsplit(x, " |_")[[1]][1:2], collapse = " "))
  tip.spp.duplicates <- which(duplicated(sppOnly.vector))
  tr.uniques$sppOnlyOrigTips <- tr.uniques$sppOnly <- drop.tip(tr.uniques$all, tip.spp.duplicates)
	tr.uniques$sppOnly$tip.label <- sppOnly.vector[!duplicated(sppOnly.vector)]
  tr.uniques$sppOnlyOrigTips$tip.label <- dat.renamed[tr.uniques$sppOnlyOrigTips$tip.label, 'tips.orig']
  return(tr.uniques)
  }


tr.uniques <- lapply(tr, makeUniques, dat)
write.tree(tr.uniques[[1]]$all, paste(trPath, trName, '.uniqueTaxa.fullNames.tre', sep = ''))
write.tree(tr.uniques[[1]]$taxaOnly, paste(trPath, trName, '.uniqueTaxa.taxonNamesOnly.tre', sep = ''))
write.tree(tr.uniques[[1]]$sppOnly, paste(trPath, trName, '.uniqueSpp.taxonNamesOnly.tre', sep = ''))
write.tree(tr.uniques[[1]]$sppOnlyOrigTips, paste(trPath, trName, '.uniqueSpp.originalTips.tre', sep = ''))
