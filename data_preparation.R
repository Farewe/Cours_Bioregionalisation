#### species data ####

a <- read.csv('./data/database/Occurences.csv', sep = ";")
levels(a$Species) <- paste0("sp", 1:nlevels(a$Species))

base2016 <- read.table("c:/r/Projets/FishBiogeography/data/BasePoisson102016_2.csv", h = T, sep =  ";")

fish.info <- unique(base2016[, c("FB_REF2015", "fresh", "brack", "salt", "AnaCat")])



fish.fresh.sp.nomig.all <- fish.all[-which(fish.all$salt=="t" & fish.all$brack == "f" & fish.all$fresh == "f"),]
fish.fresh.sp.nomig.all <- fish.fresh.sp.nomig.all[-which(fish.fresh.sp.nomig.all$AnaCat %in% 
                                                            c("anadromous", 
                                                              "amphidromous", 
                                                              "diadromous", 
                                                              "amphidromous?", 
                                                              "anadromous?", 
                                                              "catadromous",
                                                              "diadromous", 
                                                              "oceano-estuarine", 
                                                              "oceanodromous")) ,]
fish.fresh.sp.nomig.all <- droplevels(fish.fresh.sp.nomig.all)

fishdb <- a[-which(a$Species %in% fish.info$FB_REF2015[which(fish.info$salt=="t" & 
                                                               fish.info$brack == "f" &
                                                               fish.info$fresh == "f")]), ]
fishdb <- fishdb[-which(fishdb$Species %in% fish.info$FB_REF2015[which(fish.info$AnaCat %in% 
                                                                         c("anadromous", 
                                                                           "amphidromous", 
                                                                           "diadromous", 
                                                                           "amphidromous?", 
                                                                           "anadromous?", 
                                                                           "catadromous",
                                                                           "diadromous", 
                                                                           "oceano-estuarine", 
                                                                           "oceanodromous"))]), ]
fishdb1 <- fishdb[which(fishdb$Status == "native"), c(1, 2)]
fishdb2 <- fishdb[, c(1, 2)]



save(fishdb1, file = "./data/data_cours/fishdb1.RData")
save(fishdb2, file = "./data/data_cours/fishdb2.RData")





#### shapefile data ####
library(maptools)
basins <- readShapeSpatial("./data/database/basin2013_simplif")
levels(fishdb$Basin)[!(levels(fishdb$Basin) %in% levels(basins@data$BASIN))]
levels(basins@data$BASIN)[!(levels(basins@data$BASIN) %in% levels(fishdb$Basin))]
