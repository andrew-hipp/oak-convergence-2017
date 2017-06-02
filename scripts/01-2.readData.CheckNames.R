## read all data and make sure that names of species match specimen names

dat <- read.delim('../data/tables/specimens.metadata.2017-05-15.subsetted_corrected.tsv', as.is = TRUE)
dat.fields <- read.csv('../data/tables/specimens.metadata.2017-05-15.subsetted_fieldsForPaper.csv', as.is = T, row.names = 1)
dat <- dat[, row.names(dat.fields)[which(dat.fields$keepForPaper)]]
moisture.index <- read.csv('../data/tables/Index_of_Moisture_allspp_FIXED_8-5-15_ahSubsetted_edited-2017-02-21b.csv', as.is = T)
eco.all <- read.csv('../data/tables/all.eco.data.exportedFromR.2016-02-03.csv', as.is = T, row.names = 'OBJECTID')
sect.species.translate <- read.csv('../data/tables/sect.species.translate.csv', as.is = T, row.names = 1)
soils.meta <- read.delim('../data/tables/soils.metadata.2016-02-03.tsv', as.is = T)
tip.geog <- read.csv('../data/tables/spp.geog.csv', as.is = T, row.names = 1)
lf.traits <- read.delim('../data/tables/lfPhenology.2016-03-09.jcb.tsv', row.names = 1, as.is = T)

taxa <- list(wo = row.names(sect.species.translate)[grep("Quercus", sect.species.translate$Section)],
             ro = row.names(sect.species.translate)[grep("Lobatae", sect.species.translate$Section)],
             wo.mx = row.names(sect.species.translate)[intersect(grep("Quercus", sect.species.translate$Section), grep("MX", sect.species.translate$subclade))],
             ro.mx = row.names(sect.species.translate)[intersect(grep("Lobatae", sect.species.translate$Section), grep("MX", sect.species.translate$subclade))],
             ingroup = row.names(sect.species.translate)[grep("outgroup|Cerris", sect.species.translate$Section, invert = TRUE)]
             )


eco.all$Species <- paste('Quercus', eco.all$Species, sep = "_")
moisture.index$Species <- paste('Quercus', moisture.index$Species, sep = "_")
dat$cleanedSpecies <- gsub(" ", "_", dat$cleanedSpecies, fixed = TRUE)
dat$cleanedSpecies.spOnly <- label.elements(dat$cleanedSpecies, delim = "_", returnNum = 1:2, returnDelim = "_", fixed = TRUE)
dat.check <- lapply(list(moisture.index = moisture.index$Species,
                         eco.all = eco.all$Species,
                         tip.geog = row.names(tip.geog),
                         lf.traits = row.names(lf.traits),
                         sect.species.translate = row.names(sect.species.translate)),
                    function(x) {
                      out <- list(
                        missing.in.specimenData = setdiff(x, dat$cleanedSpecies.spOnly),
                        missing.in.tipData = setdiff(dat$cleanedSpecies.spOnly, x))
                      out
                    }
                  )

# leaf traits data cleanup
lf.traits$lf1 <- lf.traits$lfPhenology
lf.traits$lf1[lf.traits$lf1 == 'Evergreen'] <- 2
lf.traits$lf1[lf.traits$lf1 == 'Brevideciduous'] <- 1
lf.traits$lf1[lf.traits$lf1 == 'Deciduous'] <- 0
lf.traits$lf1[lf.traits$lf1 == 'Deciduous, Brevideciduous'] <- 0
lf.traits$lf1[lf.traits$lf1 == ''] <- NA
lf.traits$lf1 <- as.numeric(lf.traits$lf1)
lf1 <- lf.traits$lf1
names(lf1) <- row.names(lf.traits)
lf2 <- lf1
lf2[lf2 %in% 0:1] <- 0
lf2[lf2 == 2] <- 1

# moisture index means and cleanup for drought severity analyses
im.split <- split(moisture.index, moisture.index$Species)
im.means <- cbind(Im = sapply(im.split, function(x) mean(x$Im)),
                  stdev = sapply(im.split, function(x) sd(x$Im)),
				  N = sapply(im.split, function(x) dim(x)[1]),
                  sem = sapply(im.split, function(x) sd(x$Im) / sqrt(dim(x)[1])),
				  estimated = 0)
se.fail = which(is.na(im.means[, 'sem']))
im.means[se.fail, 'sem'] <- mean(im.means[, 'stdev'], na.rm = TRUE) # substitute in average sd divided by sample size (N = 1) when N = 1
im.means[se.fail, 'estimated'] <- 1
rm(im.split)
