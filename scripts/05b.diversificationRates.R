crown.taxa <- list(oaks = c('Quercus acutissima', 'Quercus rubra'),
                   am.oaks = c('Quercus alba', 'Quercus rubra'),
                   protobalanus = c('Quercus palmeri', 'Quercus tomentella'),
        				   ponticae = c('Quercus sadleriana', 'Quercus pontica'),
                   virentes = c('Quercus fusiformis', 'Quercus geminata'),
                   quercus.ss = c('Quercus lobata', 'Quercus alba'),
                   quercus.eurasia = c('Quercus robur', 'Quercus aliena'),
        				   quercus.mx = c('Quercus germana', 'Quercus turbinella'),
        				   quercus.mx.tx = c('Quercus turbinella', 'Quercus mohriana'),
        				   quercus.az = c('Quercus turbinella', 'Quercus arizonica'),
                   lobatae = c('Quercus agrifolia', 'Quercus rubra'),
                   lobatae.mx = c('Quercus crassipes', 'Quercus canbyi'))

crown.diversity <- c(oaks = 435,
                   am.oaks = 284,
                   protobalanus = 4,
        				   ponticae = 2,
        				   virentes = 7,
                   quercus.ss = 150,
                   quercus.eurasia = 25,
        				   quercus.mx = 80,
        				   quercus.mx.tx = 86,
        				   quercus.az = 8,
                   lobatae = 120,
                   lobatae.mx = 74)

names(crown.taxa) <- names(crown.diversity) <- c('Oaks (Quercus)',
                                                 'American oak clade',
                                                 'Protobalanus',
                                                 'Ponticae',
                                                 'Virentes',
                                                 'sect. Quercus s.s.',
                                                 'sect. Quercus, Eurasian lineage',
                                                 'sect. Quercus, Mexico lineage',
                                                 'sect. Quercus, Mexico + Texas lineages',
                                                 'sect. Quercus, SW U.S. lineage',
                                                 'sect. Lobatae',
                                                 'sect. Lobatae, Mexico lineage')

calib.nodes <- sapply(crown.taxa[1:4], findMRCA, tree = trs.jackknife[[1]])
crown.nodes <- sapply(crown.taxa, findMRCA, tree = trs.jackknife[[1]])

div.ages.4c.table <- as.data.frame(t(sapply(crown.nodes, function(x) quantile(sapply(trs.calib.jk.4, function(y) max(node.depth.edgelength(y)) - node.depth.edgelength(y)[x]), c(0.025, 0.5, 0.975)))))
div.ages.2c.table <- as.data.frame(t(sapply(crown.nodes, function(x) quantile(sapply(trs.calib.jk.2, function(y) max(node.depth.edgelength(y)) - node.depth.edgelength(y)[x]), c(0.025, 0.5, 0.975)))))
div.ages.2c.table$diversity <- div.ages.4c.table$diversity <- crown.diversity
div.ages.2c.table$divRate.epsilon00 <- div.ages.2c.table$divRate.epsilon09 <- div.ages.4c.table$divRate.epsilon00 <- div.ages.4c.table$divRate.epsilon09 <- NA
div.ages.2c.table$crown <- div.ages.4c.table$crown <- TRUE

## diversification rates
for(i in 1:dim(div.ages.2c.table)[1]) {
  div.ages.2c.table$divRate.epsilon00[i] <- bd.ms(time = div.ages.2c.table[i, '50%'], n = div.ages.2c.table$diversity[i], crown = div.ages.2c.table$crown[i], epsilon = 0)
  div.ages.2c.table$divRate.epsilon09[i] <- bd.ms(time = div.ages.2c.table[i, '50%'], n = div.ages.2c.table$diversity[i], crown = div.ages.2c.table$crown[i], epsilon = 0.9)
}
div.ages.2c.table.pretty <- div.ages.2c.table
div.ages.2c.table.pretty[['Crown age']] <- apply(div.ages.2c.table[, 1:3], 1, function(x) paste(round(x[2], 1), ' (', round(x[1], 1), ',', round(x[3], 1), ')', sep = ''))
div.ages.2c.table.pretty <- div.ages.2c.table.pretty[c('Crown age', "diversity","divRate.epsilon09","divRate.epsilon00")]
div.ages.2c.table.pretty$"divRate.epsilon09" <- round(div.ages.2c.table.pretty$"divRate.epsilon09", 5)
div.ages.2c.table.pretty$"divRate.epsilon00" <- round(div.ages.2c.table.pretty$"divRate.epsilon00", 5)
names(div.ages.2c.table.pretty) <- c('Crown age (0.025, 0.975 quantile)', 'Species diversity', 'Diversification rate, epsilon = 0.9', 'Diversification rate, epsilon = 0.0')
write.csv(div.ages.2c.table.pretty, 'out/diversification.ages.2c.table.pretty.2017-05-25.csv')

for(i in 1:dim(div.ages.4c.table)[1]) {
  div.ages.4c.table$divRate.epsilon00[i] <- bd.ms(time = div.ages.4c.table[i, '50%'], n = div.ages.4c.table$diversity[i], crown = div.ages.4c.table$crown[i], epsilon = 0)
  div.ages.4c.table$divRate.epsilon09[i] <- bd.ms(time = div.ages.4c.table[i, '50%'], n = div.ages.4c.table$diversity[i], crown = div.ages.4c.table$crown[i], epsilon = 0.9)
}
div.ages.4c.table.pretty <- div.ages.4c.table
div.ages.4c.table.pretty[['Crown age']] <- apply(div.ages.4c.table[, 1:3], 1, function(x) paste(round(x[2], 1), ' (', round(x[1], 1), ',', round(x[3], 1), ')', sep = ''))
div.ages.4c.table.pretty <- div.ages.4c.table.pretty[c('Crown age', "diversity","divRate.epsilon09","divRate.epsilon00")]
div.ages.4c.table.pretty$"divRate.epsilon09" <- round(div.ages.4c.table.pretty$"divRate.epsilon09", 5)
div.ages.4c.table.pretty$"divRate.epsilon00" <- round(div.ages.4c.table.pretty$"divRate.epsilon00", 5)
names(div.ages.4c.table.pretty) <- c('Crown age (0.025, 0.975 quantile)', 'Species diversity', 'Diversification rate, epsilon = 0.9', 'Diversification rate, epsilon = 0.0')
write.csv(div.ages.4c.table.pretty, 'out/diversification.ages.4c.table.pretty.2017-05-25.csv')
