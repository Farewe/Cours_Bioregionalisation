# 1. 
library(biogeonetworks)
library(rgdal)
library(plyr)
library(RColorBrewer)
# Shapefiles
basins <- readOGR("./data/data_cours/basin2013_simplif.shp")

# Occurrence databases
load("./data/data_cours/fishdb1.RData")
load("./data/data_cours/fishdb2.RData")


# 2.
writePajek(fishdb1, 
           site.field = "Basin",
           species.field = "Species",
           filename = "./data/fishdb1.net")


# 3.
system("./infomap --undirected --tree --map ./data/fishdb1.net ./data/")


# 4.
fish.net1 <- readInfomapTree("./data/fishdb1.tree")

# 5. 
fish.sites1 <- getSiteTable(fishdb1, site.field = "Basin", network = fish.net1)
fish.species1 <- getSpeciesTable(fishdb1, species.field = "Species", network = fish.net1)

head(fish.sites1)
count(fish.sites1$lvl1)
count(fish.species1$lvl1)

# 6. 
fish.net1 <- attributeColors(fish.net1, nb.max.colors = 9, db = fishdb1)
fish.sites1 <- attributeColors(fish.sites1, nb.max.colors = 9, db = fishdb1)


# 7.
basins@data$color1 <- fish.net1$color[match(basins@data$BASIN, fish.net1$Name)]
plot(basins, col = basins@data$color1)


# 8. 
writeGDF(fishdb1, fish.net1, site.field = "Basin", species.field = "Species", color.field = "color",
         filename = "./data/fishnet1.gdf")

# 10. 
writePajek(fishdb2, 
           site.field = "Basin",
           species.field = "Species",
           filename = "./data/fishdb2.net")

system("./infomap --undirected --tree --map ./data/fishdb2.net ./data/")

fish.net2 <- readInfomapTree("./data/fishdb2.tree")


fish.sites2 <- getSiteTable(fishdb2, site.field = "Basin", network = fish.net2)
fish.species2 <- getSpeciesTable(fishdb2, species.field = "Species", network = fish.net2)

head(fish.sites2)
count(fish.sites2$lvl1)
count(fish.species2$lvl1)

fish.net2 <- attributeColors(fish.net2, nb.max.colors = 6, db = fishdb2)
fish.sites2 <- attributeColors(fish.sites2, nb.max.colors = 6, db = fishdb2)

basins@data$color2 <- fish.net2$color[match(basins@data$BASIN, fish.net2$Name)]

op <- par(mfrow = c(2, 1), mar = c(1.1, 1.1, 1.1, 1.1))
plot(basins, col = basins@data$color1)
plot(basins, col = basins@data$color2)
par(op)

writeGDF(fishdb2, fish.net2, site.field = "Basin", species.field = "Species", color.field = "color",
         filename = "./data/fishnet2.gdf")

# 11.
fish.net1 <- participationCoefficient(fish.net1, 
                                      db = fishdb1, site.field = "Basin",
                                      species.field = "Species",
                                      lvl = "lvl1")
fish.net2 <- participationCoefficient(fish.net2, 
                                      db = fishdb2, site.field = "Basin",
                                      species.field = "Species",
                                      lvl = "lvl1")
rbPal <- colorRampPalette(c('blue','red'))
fish.net1$pcCol <- rbPal(10)[as.numeric(cut(fish.net1$participation.coef, breaks = 10))]
basins@data$pcCol1 <- fish.net1$pcCol[match(basins@data$BASIN, fish.net1$Name)]
fish.net2$pcCol <- rbPal(10)[as.numeric(cut(fish.net2$participation.coef, breaks = 10))]
basins@data$pcCol2 <- fish.net2$pcCol[match(basins@data$BASIN, fish.net2$Name)]

op <- par(mfrow = c(2, 1), mar = c(1.1, 1.1, 1.1, 1.1))
plot(basins, col = basins@data$pcCol1)
plot(basins, col = basins@data$pcCol2)
par(op)


