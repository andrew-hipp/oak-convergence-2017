## Originally summer 2016
## updated 2017-01-04 to add bootstraps
## completely rewritten 2017-04-21 for so many reasons

require("rgdal") # requires sp, will use proj.4 if installed
require("maptools")
require("ggplot2")
require("plyr")
require('dismo')
require('raster')
require('grid')
require(geoscale)
data(timescales)

## new plot starts here

if(!exists('path')) path <- 'out'
pdf(paste(path, '/FIG02.timetree.newVersion.pdf', sep = ''), 7, 9.5)
par(mar = c(2.4,0,0,0) + 0.1)

regionColors = c(C = 'deepskyblue', E = 'maroon', M = 'orange', U = 'black')

cladesOfInterest = list(
   root = c('Lithocarpus_hancei', 'Quercus_arizonica'),
   newWorldOaks = c('Quercus_kelloggii', 'Quercus_arizonica'),
   Cerris=c('Quercus_baronii', 'Quercus_trojana'),
   Lobatae.core = c('Quercus_kelloggii', 'Quercus_gentryi'),
   Lobatae.Mexico = c('Quercus_canbyi', 'Quercus_hypoleucoides'),
   Protobalanus = c('Quercus_palmeri', 'Quercus_tomentella'),
   Ponticae = c('Quercus_sadleriana', 'Quercus_pontica'),
   Virentes = c('Quercus_fusiformis', 'Quercus_minima'),
   Quercus.core = c('Quercus_lobata', 'Quercus_alba'),
   Quercus.Eurasia = c('Quercus_mongolica', 'Quercus_robur'),
   Quercus.Mexico = c('Quercus_potosina', 'Quercus_germana'),
   Quercus.TX = c('Quercus_pungens', 'Quercus_polymorpha'),
   Quercus.AZ = c('Quercus_turbinella', 'Quercus_arizonica')
  )

  overSections = list(
    Cerris=c('Quercus_baronii', 'Quercus_trojana'),
    Lobatae = c('Quercus_kelloggii', 'Quercus_hypoleucoides'),
    Protobalanus = c('Quercus_palmeri', 'Quercus_tomentella'),
    Ponticae = c('Quercus_sadleriana', 'Quercus_pontica'),
    Virentes = c('Quercus_fusiformis', 'Quercus_minima'),
    "Texas / N. Mexico" = c('Quercus_mohriana', 'Quercus_vaseyana')
    )

  underSections = list(
    "section Quercus" = c('Quercus_lobata', 'Quercus_arizonica'),
    "American oaks" = c('Quercus_lobata', 'Quercus_kelloggii'),
    "Arizona / N. Mexico" = c('Quercus_laeta', 'Quercus_turbinella')
    )


  # get tree, plot with time scale
  t1 <- subset(timescales$ICS2013, timescales$ICS2013[, 'Type'] == 'Epoch')
  tr.beast <- tr.spp.4c.discreteClock.beast
  tr.temp <- tr.spp.4c.discreteClock

  a=plot(tr.temp, edge.color = 0, tip.color = 0, x.lim = c(-7, 96.61137))
  #abline(v = 66.34929-t1[2:7, 'Start'], lty = 'dashed', col = 'gray')
  segments(x0 = 66.34929-t1[2:7, 'Start'], y0 = 0, y1 = length(tr.temp$tip.label) + 1, lty = 'dashed', col = 'gray')
  rect(66.34929-t1[2:7, 'Start'], par()$usr[3]+2, 66.34929-t1[2:7, 'End'], par()$usr[3] + 4, col = c('white', 'gray80'), border = 'black')
  text(66.34929-t1[2:7, 'Start'], par()$usr[3]+1, labels = round(t1[2:7, 'Start'], 1), cex = 0.5)
  text(66.34929-t1[2:7, 'Midpoint'], par()$usr[3] - 1.5, labels = t1[2:7, 'Name'], srt = -60, adj = 0, xpd = TRUE, cex = 0.5)
  plot.phylo.upon(tr.temp, cex = 0.4)

  # add biogeography coding
  t2 <- tip.geog[tr.temp$tip.label,c('C', 'E', 'M', 'U')]
  for(i in seq(dim(t2)[[1]])){
    for(j in seq(dim(t2)[2])) {
      points(j + 83, i, pch = 22, bg = ifelse(t2[i,j] == 1, regionColors[j], "white"), lwd = 0.5)
    }}
  text(84:(83+dim(t2)[2]), rep(length(tr.temp$tip.label) + 2, dim(t2)[2]), dimnames(t2)[[2]], cex = 0.5, xpd = TRUE)

  # add leaf phenology
  t3 <- lf.traits[tr.temp$tip.label, 'lfPhenology']
  t3[t3 == "Deciduous, Brevideciduous"] <- 'Deciduous'
  t3[is.na(t3)] <- ''
  t3[t3 == ''] <- 'notCoded'
  t3.colors <- sapply(t3, function(x) switch(x,
                                             Deciduous = 'white',
                                             Brevideciduous = 'lightgreen',
                                             Evergreen = 'darkgreen',
                                             notCoded = NA))
  points(rep(81, length(t3)), 1:length(t3), pch = ifelse(t3 == 'notCoded', NA, 22), bg = t3.colors, lwd = 0.4)
  text(81, length(tr.temp$tip.label) + 2.2, "Leaf\nhabit", cex = 0.4)

  # add node bars
  tr.mrca <- mrca(tr.temp)
  nodesToBar <- sapply(cladesOfInterest, function(x) tr.mrca[x[1], x[2]])
  par(lend = 2)
  HPDbars(tr.beast, nodes = nodesToBar, col = 'black', lty = 'solid', lwd = 4.5)
  HPDbars(tr.beast, nodes = nodesToBar, col = 'gray', lty = 'solid', lwd = 4)

  ## section / clade names
  lastP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  nodesToSection <- sapply(overSections, function(x) tr.mrca[x[1], x[2]])
  for(i in 1:length(nodesToSection)) text(lastP$xx[nodesToSection[i]],
                                        lastP$yy[nodesToSection[i]],
                                        names(nodesToSection)[i],
                                        adj = c(1.1,-0.5),
                                        offset = c(0,1.2),
                                        cex = 0.6,
                                        font = 2)

  nodesToSection <- sapply(underSections, function(x) tr.mrca[x[1], x[2]])
  for(i in 1:length(nodesToSection)) text(lastP$xx[nodesToSection[i]],
                                          lastP$yy[nodesToSection[i]],
                                          names(nodesToSection)[i],
                                          adj = c(1.1,1.7),
                                          offset = c(0,-1.2),
                                          cex = 0.6,
                                          font = 2)


  ## leaf habit legend
  legend(x=-6.5, y=46,
         cex = 0.5, pch = 22, pt.cex = 1.5,
         pt.bg = c('darkgreen', 'lightgreen', 'white'),
         legend = c("Evergreen", "Brevideciduous", "Deciduous"),
         box.lwd = 0.5, box.col = 'white', bg = 'white',
         title = "")
  legend(x = -6.5, y = 46.5, legend = c('','', ''),
    title = 'Leaf habit', title.adj = 0.15,
    bty = 'n', cex = 0.6)

  ## map legend
  vp1 <- viewport(x = 0.03, y = 0.09,
                  width = 0.33, height = 0.20,
                  just = c('left', 'bottom'),
                  gp = gpar(bty = 'o'))
  b = ggplot(mi.subsampled.bySp, aes(x=lon, y=lat))
  ourArea.noCanada <- map_data('world', c('usa','mexico',
                                 'costa rica', 'honduras', 'guatemala', 'belize', 'el salvador', 'nicaragua', 'panama'
                                 ))
  ourArea.states <- map_data('state')
  b <- b + geom_polygon(data = regions.df, aes(long, lat, group = group, fill = oakAreas))
  b <- b + scale_fill_manual('Biogeographic\nregions',
                             labels = c('C', 'E', 'M'),
                             values = as.character(regionColors[1:3]))
  b = b + geom_map(data = ourArea.noCanada, map = ourArea.noCanada, aes(x = long, y = lat, map_id = region), colour = 'black', fill = NA, size = 0.2)
  b = b + geom_map(data = ourArea.states, map = ourArea.states, aes(x = long, y = lat, map_id = region), colour = 'black', fill = NA, size = 0.1)
  b = b + xlim(-125, -57) + ylim(8, 50)
  # for(i in seq(length(reg.list))) b <- b + geom_path(data = regions.df[regions.df$bioregio %in% reg.list[[i]], ], aes(long, lat, group = group), colour = reg.colors.under[i], size = 0.1)
  #b <- b + geom_path(data = regions.df, aes(long, lat, group = group, colour = oakAreas), size = 0.3)
  # b <- b + scale_colour_manual("Biogeographic\nregions", values=reg.colors.under)

  b = b + theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank())
  b <- b + theme(legend.position = c(0.85, 0.3),
                 legend.background = element_rect(fill=NA),
                 legend.title = element_text(size = 0),
                 legend.text = element_text(size = 5),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_rect(colour = 'gray', fill = NA, size = 0.5),
                panel.background = element_blank(),
                plot.title = element_text(size = 7)
                 )
  b <- b + ggtitle("Biogeographic regions")
  print(b, vp = vp1)

#

dev.off()
