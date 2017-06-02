## clean up data and write bioregions file
## also reads and preps ecology data

clean.eco <- TRUE
make.maps <- TRUE
summarize.sampling <- TRUE
plot.latitude <- TRUE
plot.by.rank <- FALSE
mapPath = 'out/species.maps.ecoPoints.postEdit.v4/'

if(length(dir(mapPath)) == 0) dir.create(mapPath)

if(clean.eco) {
  moisture.index$thinField <- paste(round(moisture.index$longitude, 2),
                                    round(moisture.index$latitude, 2),
                                    moisture.index$Species,
                                    sep = '_')
  moisture.index.subsetted <- moisture.index[which(moisture.index$use), ]
  moisture.index.thinned <- moisture.index.subsetted[!duplicated(moisture.index.subsetted$thinField), ]
  eco.oaks <- cbind(moisture.index.thinned,
                    eco.all[as.character(moisture.index.thinned$OBJECTID), -which(names(eco.all) %in% names(moisture.index.thinned))])
  eco.oaks <- eco.oaks[c("Species",
                         grep('bio|latitude|longitude', names(eco.oaks), value = TRUE),
                         soils.meta$FIELD[which(soils.meta$utilize == 1)])]
  eco.oaks$T_TEXTURE[which(eco.oaks$T_TEXTURE == 0)] <- NA
  eco.split <- split(eco.oaks[-1], eco.oaks$Species)
  #names(eco.split) <- paste('Quercus', names(eco.split), sep = '_')
  eco.means <- t(sapply(eco.split, function(x) apply(x, 2, mean, na.rm = T)))
  eco.sd <- t(sapply(eco.split, function(x) apply(x, 2, sd, na.rm = T)))
  eco.N <- sapply(eco.split, function(x) dim(x)[1])
  # eco.means[which(eco.N == 1), grep('bio|latitude|longitude', names(eco.means))] <- do.call(rbind, eco.split[which(eco.N == 1)])[, grep('bio|latitude|longitude', names(eco.means))]
  eco.spp <- intersect(tr.spp.4c.discreteClock.noOG$tip.label, row.names(eco.means))
  tr.eco <- drop.tip(tr.spp.4c.discreteClock.noOG, which(!tr.spp.4c.discreteClock.noOG$tip.label %in% eco.spp))
  eco.means <- eco.means[eco.spp, ]
  biosort <- order(eco.means[, 'latitude'], decreasing = FALSE)
  # biosort[which(eco.N < 5)] <- NULL
  eco.means <- eco.means[biosort, ]
  eco.sd <- eco.sd[biosort, ]
  eco.N <- eco.N[biosort]

  eco.means <- list(soil = eco.means[ , which(colnames(eco.means) %in% soils.meta$FIELD[soils.meta$utilize == 1])],
                  clim = eco.means[ , grep('bio', colnames(eco.means), value = TRUE)],
                  all = eco.means)

  write.csv(moisture.index.thinned[ , c('Species', 'latitude', 'longitude')], paste('bioregions.table', Sys.Date(), 'csv', sep = '.'))

}

if(make.maps) {
  for(i in unique(moisture.index.thinned$Species)){
    pdf(paste(mapPath, i, '.pdf', sep = ''))
    map('world', c('usa', 'mexico', 'costa rica', 'honduras', 'guatemala', 'belize', 'el salvador', 'nicaragua', 'panama', 'colombia'), xlim = c(-130, -60))
    points(moisture.index.thinned$longitude[moisture.index.thinned$Species == i][which(moisture.index.thinned$use[moisture.index.thinned$Species == i])],
           moisture.index.thinned$latitude[moisture.index.thinned$Species == i][which(moisture.index.thinned$use[moisture.index.thinned$Species == i])],
           pch = 19, col = 'red')
    abline(h = seq(from = 0, to = 120, by = 10), lty = 'dashed')
    abline(v = seq(from = -120, to = -70, by = 10), lty = 'dashed')
    text(x = c(rep(-125, 13), seq(from = -120, to = -70, by = 10))-2,
         y = c(seq(from = 0, to = 120, by = 10), rep(10, 6))+1,
         label = c(seq(from = 0, to = 120, by = 10), seq(from = -120, to = -70, by = 10)),
         cex = 0.5)
    text(x = c(-89.8, -93.2), y = c(40.0, 45.0), label = c('STRHL', 'UMN'), cex = 0.5)
    dev.off()
  }
}

if(summarize.sampling) {
  message(paste("total spp in specimen table :", length(unique(moisture.index.thinned$Species))))
  message(paste("total specimens :", dim(moisture.index.thinned)[1]))
  mis.s <- sapply(split(moisture.index.thinned, moisture.index.thinned$Species), function(x) dim(x)[1])
  message(paste("specimens / sp :", mean(mis.s), "+/-", sd(mis.s), "(sd)"))
  message("specimen range --")
  print(range(mis.s))
  message(paste('spp with 10 or fewer records :', sum(mis.s < 11)))
  message(paste('spp with at least 50 records :', sum(mis.s > 49)))
  message(paste('spp with at least 100 records :', sum(mis.s > 99)))
  print(tail(sort(mis.s)))
}

if(plot.latitude) {
  pdf('out/eco.means.plot.v5.pdf')
  if(plot.by.rank) {x <- seq(dim(eco.means$all)[1])}
    else x <- eco.means$all[, 'latitude']
  plot(x, eco.means$all[, 'bio11'] / 10, type = 'l', col = 'blue', lwd =2,
       ylim = range(eco.means$all[, c('bio11', 'bio10', 'bio6')] / 10, na.rm = T),
       xlab = paste("Mean latitude for each of", dim(eco.means$all)[1], "species"),
       ylab = "Mean temperature for each species (degrees Celsius)")
  lines(x, eco.means$all[, 'bio10'] / 10, type = 'l', col = 'red', lwd = 2)
  lines(x, eco.means$all[, 'bio6'] / 10, type = 'l', col = 'blue', lwd = 1, lty = 'dashed')
  legend(10, -10, legend = c('BIO10 - Mean temp, warmest quarter',
                          'BIO11 - Mean temp, coldest quarter',
                          'BIO6  - Min temp, coldest month'),
          cex = 0.7, lwd = c(2, 2, 1), lty = c('solid', 'solid', 'dashed'),
          col = c('red', 'blue', 'blue'),
          bty = 'n')
  dev.off()
}
