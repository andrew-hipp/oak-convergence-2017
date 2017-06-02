source('../scripts/00.taxa.R')
eco.wo <- c("Quercus_engelmannii", "Quercus_oblongifolia", "Quercus_laeta",
"Quercus_arizonica", "Quercus_ajoensis", "Quercus_turbinella",
"Quercus_toumeyi", "Quercus_grisea", "Quercus_glaucoides", "Quercus_peduncularis",
"Quercus_deserticola", "Quercus_potosina", "Quercus_obtusata",
"Quercus_diversifolia", "Quercus_purulhana", "Quercus_resinosa",
"Quercus_liebmanii", "Quercus_glaucescens", "Quercus_martinezii",
"Quercus_rugosa", "Quercus_greggii", "Quercus_glabrescens", "Quercus_insignis",
"Quercus_corrugata", "Quercus_lancifolia", "Quercus_germana",
"Quercus_pungens", "Quercus_vaseyana", "Quercus_laceyi", "Quercus_mohriana",
"Quercus_austrina", "Quercus_chapmanii", "Quercus_havardii",
"Quercus_margarettae", "Quercus_stellata", "Quercus_similis",
"Quercus_boyntonii", "Quercus_oglethorpensis", "Quercus_sinuata",
"Quercus_macrocarpa", "Quercus_lyrata", "Quercus_bicolor", "Quercus_prinoides",
"Quercus_muehlenbergii", "Quercus_alba", "Quercus_montana", "Quercus_douglasii",
"Quercus_john-tuckeri", "Quercus_dumosa", "Quercus_pacifica",
"Quercus_durata", "Quercus_berberidifolia", "Quercus_garryana",
"Quercus_lobata")

eco.ro <- c("Quercus_hypoleucoides", "Quercus_scytophylla", "Quercus_candicans",
"Quercus_crassifolia", "Quercus_conzattii", "Quercus_polymorpha",
"Quercus_conspersa", "Quercus_emoryi", "Quercus_peninsularis",
"Quercus_eduardii", "Quercus_gentryi", "Quercus_crassipes", "Quercus_castanea",
"Quercus_laurina", "Quercus_pinnativenulosa", "Quercus_crispifolia",
"Quercus_affinis", "Quercus_acutifolia", "Quercus_mexicana",
"Quercus_cortesii", "Quercus_costaricensis", "Quercus_eugeniifolia",
"Quercus_sapotifolia", "Quercus_segoviensis", "Quercus_elliptica",
"Quercus_iltisii", "Quercus_uxoris", "Quercus_planipocula", "Quercus_canbyi",
"Quercus_myrtifolia", "Quercus_inopina", "Quercus_incana", "Quercus_hemisphaerica",
"Quercus_laevis", "Quercus_marilandica", "Quercus_nigra", "Quercus_laurifolia",
"Quercus_elliottii", "Quercus_phellos", "Quercus_imbricaria",
"Quercus_georgiana", "Quercus_arkansana", "Quercus_falcata",
"Quercus_pagoda", "Quercus_ilicifolia", "Quercus_shumardii",
"Quercus_acerifolia", "Quercus_buckleyi", "Quercus_rubra", "Quercus_coccinea",
"Quercus_ellipsoidalis", "Quercus_velutina", "Quercus_palustris",
"Quercus_texana", "Quercus_wislizeni", "Quercus_parvula", "Quercus_agrifolia",
"Quercus_kelloggii")

pchSect <- c(21, 25, 22, 22, 23, 24)
colorSect <- c('red', 'white', 'green', 'blue', 'black', 'yellow')
names(colorSect) <- names(pchSect) <- c('Lobatae', 'Protobalanus', 'Quercus - American', 'Quercus - Eurasian', 'Sadlerianae', 'Virentes')

colorRegion <- c('blue', 'black', 'red')
names(colorRegion) <- c("ENA", "MX.SW", "CA")
