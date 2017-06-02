## make a decent map of our species records

## may require rgeos and gpclib if gpclibPermitStatus == FALSE

dat.renamed.cleaned <- read.csv('../data/tables/dat.renamed.cleaned.csv', row.names = 1, as.is = T)
dat.renamed.cleaned$Section <- as.factor(dat.renamed.cleaned$Section)
makeData = F
resolution = 50
from = 'white'
to = 'black'

if(makeData) {
  regions <- readOGR(dsn = '../data/bioregions.out.v2/bioregions.table.2017-02-22')
  regions@data$id = rownames(regions@data)
  regions.points = fortify(regions, region="id")
  regions.df = join(regions.points, regions@data, by="id")
  regions.df$oakAreas = ""
  reg.list <- list("ENA" = c(1,2,6, 9),
                   "CA-FP" = c(3,11),
                   "MX" = c(4,5,7,8,10))
  for(i in names(reg.list)) regions.df$oakAreas[regions.df$bioregio %in% reg.list[[i]]] <- i
  }

if(makeData) {
  mi.subsampled.bySp <- NULL
  ourArea <- map_data('world', c('usa', 'canada', 'mexico',
                                 'costa rica', 'honduras', 'guatemala', 'belize', 'el salvador', 'nicaragua', 'panama'
                                 ))
  latlon <- data.frame(lon = moisture.index.thinned$longitude, lat = moisture.index.thinned$latitude)
  rr <- raster(extent(range(latlon$lon), range(latlon$lat)), nrows = resolution, ncols = resolution)
  for(i in unique(moisture.index.thinned$Species)) {
    rows <- which(moisture.index.thinned$Species == i)
    message(paste('doing',i, 'with', length(rows), 'specimens'))
    temp.ll <- cbind(lon = moisture.index.thinned$longitude[rows], lat = moisture.index.thinned$latitude[rows])
    if(length(rows) == 1) mi.subsampled.bySp <- rbind(mi.subsampled.bySp, temp.ll)
    else mi.subsampled.bySp <- rbind(mi.subsampled.bySp, gridSample(temp.ll, rr))
  }
  mi.subsampled.bySp <- as.data.frame(mi.subsampled.bySp)
}


pdf('out/mi.by.specimen.pdf')
a = ggplot(latlon, aes(x=lon, y=lat)) + stat_bin2d(bins = resolution) + scale_fill_gradient(low = from, high = to)
a = a + geom_map(data = ourArea, map = ourArea, aes(x = long, y = lat, map_id = region), colour = 'black', fill = NA)
a = a + xlim(-125, -63) + ylim(8, 55)
a = a + xlab("Longitude") + ylab('Latitude')
print(a)
dev.off()
# mi.subsampled.bySpecimen <- gridSample(latlon, raster(extent(range(latlon$lon), range(latlon$lat)), nrows = 100, ncols = 100))
# ggplot(mi.subsampled.bySpecimen, aes(x=lon, y=lat)) + stat_bin2d(bins = 99) + scale_fill_gradient(low = "lightblue", high ="purple4")

#pdf('out/mi.by.species.v18.alpha85-No-viewport.v10.pdf', 9, 8)
if(!exists('path')) path <- 'out'
pdf(paste(path, '/FIG01.samplingMap.pdf', sep = ''), 7, 5)
par(mar = rep(0,4))
reg.colors.under = c('red', 'green', 'blue')
reg.colors.solids = c('gray65', 'gray75', 'gray85')
reg.colors.over = c('gray70', 'maroon2', 'lightgreen')
q.colors = c(Quercus = 'black', Lobatae = 'red', Protobalanus = 'yellow')
q.shapes = c(Quercus = 17, Lobatae = 19, Protobalanus = 18)
q.sizes = c(Quercus = 2, Protobalanus = 4, Lobatae = 2)

a = ggplot(mi.subsampled.bySp, aes(x=lon, y=lat))
a = a + stat_bin2d(bins = resolution, alpha = 0.85)
a = a + scale_fill_gradient(low = from, high = to)
a = a + geom_map(data = ourArea, map = ourArea, aes(x = long, y = lat, map_id = region), colour = 'black', fill = NA)
a = a + xlim(-125, -63) + ylim(8, 49)
a = a + xlab("Longitude") + ylab('Latitude')
a = a + labs(fill = "Oak species\nper grid cell")
dat.renamed.cleaned$Long <- as.numeric(dat.renamed.cleaned$Long)
dat.renamed.cleaned$Lat <- as.numeric(dat.renamed.cleaned$Lat)
#a = a + geom_path(data = regions.df, aes(long, lat, group = group), colour = "gray", size = 0.5)
# for(i in seq(length(reg.list))) a <- a + geom_path(data = regions.df[regions.df$bioregio %in% reg.list[[i]], ], aes(long, lat, group = group), colour = reg.colors.over[i], size = 0.35)
a = a + geom_point(data = dat.renamed.cleaned[dat.renamed.cleaned$Section %in% c('Quercus', 'Lobatae', 'Protobalanus'), ],
                  mapping=aes(x=Long, y=Lat, shape=factor(Section), color = factor(Section), size = factor(Section)),
                  stroke = 0.5)
#for(i in c('Lobatae', 'Quercus', 'Protobalanus')) {
#  a = a + geom_point(data = dat.renamed.cleaned[dat.renamed.cleaned$Section == i, ],
#                     mapping=aes(x=Long, y=Lat), shape = q.shapes[i], size = q.sizes[i],
#                     colour = 'black', fill = q.colors[i], stroke = 0.5)
#  }
a = a + scale_color_manual(values=q.colors, name = "DNA Samples\nby clade")
a = a + scale_shape_manual(values=q.shapes, name = "DNA Samples\nby clade")
a = a + scale_size_manual(values=q.sizes, name = "DNA Samples\nby clade")
#a = a + scale_shape_discrete(solid = TRUE)
# a = a + geom_point(data = data.frame(lon = -90.766667, lat = 79.433333), mapping=aes(x = lon, y = lat), shape = 10, size = 5, colour = 'black')
a = a + theme(legend.position = c(0.87, 0.3),
              legend.background = element_rect(fill=NA),
              legend.key.size = unit(0.4, 'cm'))
a = a + coord_equal()
print(a)
dev.off()
