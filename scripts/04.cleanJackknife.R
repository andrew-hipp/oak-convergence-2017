## adapted from cleanTree to deal with jackknife trees
## updated 2016-01-04 for new files

ogRenamer <- c('>DM46', '>DM71', '>DM74')
trPath <- ('../data/phylogeny/rellSearch/jackknife.m20/out')

## read data
trs.jackknife <- lapply(dir(trPath, full = T, patt = 'result'), read.tree)
class(trs.jackknife) <- 'multiPhylo'
# dat <- read.delim('../data/tables/plate.names.2016.11.30.tsv', as.is = T)

## do some renaming and export trees
for(trNum in 1:length(trs.jackknife)) {
  trs.jackknife[[trNum]] <- ladderize(root(trs.jackknife[[trNum]], ogRenamer))
  # trs.jackknife[[trNum]] <- drop.tip(trs.jackknife[[trNum]], ogRenamer)
	tips.orig <- trs.jackknife[[trNum]]$tip.label

  trs.jackknife[[trNum]]$tip.label <- gsub(">", "", trs.jackknife[[trNum]]$tip.label)
	trs.jackknife[[trNum]]$tip.label <- gsub(".barcodeStripped", "", trs.jackknife[[trNum]]$tip.label)
	trs.jackknife[[trNum]]$tip.label <- gsub(".nameFixed", "", trs.jackknife[[trNum]]$tip.label)
	trs.jackknife[[trNum]]$tip.label <- gsub(".techRep", "", trs.jackknife[[trNum]]$tip.label)

  tips.clean <- trs.jackknife[[trNum]]$tip.label

	tipsIndex <- match(trs.jackknife[[trNum]]$tip.label, dat$tip.label)
  tip.trans <- structure(paste(dat$cleanedSpecies[tipsIndex], trs.jackknife[[trNum]]$tip.label, dat$geoPol[tipsIndex], sep = ' | '), names = trs.jackknife[[trNum]]$tip.label)
	dat.renamed <- dat[tipsIndex, ]
	row.names(dat.renamed) <- as.character(tip.trans)
	trs.jackknife[[trNum]]$tip.label <- as.character(tip.trans)
  dat.renamed$tips.orig <- tips.orig

  sppOnly.vector <- sapply(trs.jackknife[[trNum]]$tip.label, function(x) paste(strsplit(x, " |_", fixed = function(x) {})[[1]][1:2], collapse = " "))
  trs.jackknife[[trNum]]$tip.label <- sppOnly.vector
  }
