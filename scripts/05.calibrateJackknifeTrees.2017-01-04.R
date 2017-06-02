## analyses post-tree
## 2016-05-10 ff


crown.taxa <- list(oaks = c('Quercus acutissima', 'Quercus rubra'),
                   am.oaks = c('Quercus alba', 'Quercus rubra'),
                   lobatae = c('Quercus agrifolia', 'Quercus rubra'),
                   quercus.ss = c('Quercus lobata', 'Quercus alba'),
				   ponticae = c('Quercus sadleriana', 'Quercus pontica'),
				   quercus.eurasia = c('Quercus robur', 'Quercus aliena'),
				   quercus.mx = c('Quercus germana', 'Quercus turbinella'),
				   quercus.mx.tx = c('Quercus turbinella', 'Quercus mohriana'),
				   quercus.az = c('Quercus turbinella', 'Quercus arizonica'),
				   lobatae.mx = c('Quercus crassipes', 'Quercus canbyi'))

crown.diversity <- c(oaks = 425,
                   am.oaks = 284,
                   lobatae = 120,
                   quercus.ss = 150,
				   ponticae = 2,
				   quercus.eurasia = 25,
				   quercus.mx = 80,
				   quercus.mx.tx = 86,
				   quercus.az = 8,
				   lobatae.mx = 74)

calib.nodes <- sapply(crown.taxa[1:4], findMRCA, tree = trs.jackknife[[1]])
crown.nodes <- sapply(crown.taxa, findMRCA, tree = trs.jackknife[[1]])

roll.the.dice <- function(x = NA) {
  out <- c(oaks = rexp(1, rate = 0.3) + 48,
              am.oaks = rnorm(1, mean = 45, sd = 3/1.96),
              lobatae = rexp(1, rate = 0.3) + 31,
              quercus.ss = rexp(1, rate = 0.3) + 31)
  if(!(out[1] > out[2] & out[2] > out[3] & out[2] > out[4])) {
    print('...recursing...')
    out <- roll.the.dice()
  }
  return(out)
}

calibs.500 <- sapply(1:500, roll.the.dice)

## beware! this will change if the tree is changed... set here now for convenience. Updated 2016-05-10

calibs.tree.index <- cbind(tree = rep(1:100, 5),
                           calib = 1:500)

calibrate.a.tree <- function(x, tr = trs.jackknife, calib = calibs.500, index = calibs.tree.index, which.use = length(calib.nodes)) {
  tr.calib <- structure(list(node = as.integer(calib.nodes)[which.use],
                                    age.min = as.numeric(calib[which.use, x]),
                                    age.max = as.numeric(calib[which.use, x]),
                                    soft.bounds = rep(FALSE, length(which.use))),
                                    row.names = which.use,
                                    .Names = c("node", "age.min", "age.max", "soft.bounds"),
                                    class = "data.frame")

  tr.out <- chronos(tr[[index[x, 'tree']]], model = 'discrete', calibration = tr.calib)
  tr.out
}

trs.calib.jk.4 <- mclapply(1:500, calibrate.a.tree, which.use = 1:4, mc.cores = 14)
trs.calib.jk.2 <- mclapply(1:500, calibrate.a.tree, which.use = 1:2, mc.cores = 14)

class(trs.calib.jk.4) <- class(trs.calib.jk.2) <- 'multiPhylo'
write.tree(trs.calib.jk.4, 'out/trs.calib.jackknife.4.tre')
write.tree(trs.calib.jk.2, 'out/trs.calib.jackknife.2.tre')

save(trs.calib.jk.4, trs.calib.jk.2, calibs.500, trs.jackknife, calibs.tree.index, calib.nodes, crown.taxa, calibrate.a.tree, file = 'out/trs.calib.jackknife.Rdata')

message("YOU NEED TO RUN THESE TREE FILES THROUGH TREE ANNOTATOR TO GET A TREE YOU CAN DRAW AND DO THINGS WITH")
message("After making them, read annoated files back in as tr.spp.4c.discreteClock and tr.spp.2c.discreteClock, /nusing read.beast from the phyloch package")
