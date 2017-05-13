# 1. 
source("./scripts/network_functions.R")

# Shapefiles
basins <- readShapeSpatial("./data/data_cours/basin2013_simplif")

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
fish.net1 <- read.infomap.tree("./data/fishdb1.tree")

# 5. 
fish.sites1 <- getSiteTable(fishdb1, site.field = "Basin", network = fish.net1)
fish.species1 <- getSpeciesTable(fishdb1, species.field = "Species", network = fish.net1)

head(fish.sites1)
summary(fish.sites1$lvl1)
summary(fish.species1$lvl1)

# 6. 
fish.net1 <- attributeColors(fish.net1, nb.max.colors = 9)
fish.sites1 <- attributeColors(fish.sites1, nb.max.colors = 9)
fish.species1 <- attributeColors(fish.species1, nb.max.colors = 9)

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

fish.net2 <- read.infomap.tree("./data/fishdb2.tree")


fish.sites2 <- getSiteTable(fishdb2, site.field = "Basin", network = fish.net2)
fish.species2 <- getSpeciesTable(fishdb2, species.field = "Species", network = fish.net2)

head(fish.sites2)
summary(fish.sites2$lvl1)
summary(fish.species2$lvl1)

fish.net2 <- attributeColors(fish.net2, nb.max.colors = 6)
fish.sites2 <- attributeColors(fish.sites2, nb.max.colors = 6)
fish.species2 <- attributeColors(fish.species2, nb.max.colors = 6)

basins@data$color2 <- fish.net2$color[match(basins@data$BASIN, fish.net2$Name)]

op <- par(mfrow = c(2, 1), mar = c(1.1, 1.1, 1.1, 1.1))
plot(basins, col = basins@data$color1)
plot(basins, col = basins@data$color2)
par(op)

writeGDF(fishdb2, fish.net2, site.field = "Basin", species.field = "Species", color.field = "color",
         filename = "./data/fishnet2.gdf")

# 11.
plotRadial(fish.sites1, levels = grep("lvl", colnames(fish.sites1)), leaf.names = "Name")
library(plyr)
fish.net1 <- partic.coef(fish.net1, 
                         db = fishdb1, site.field = "Basin",
                         species.field = "Species",
                         cluster.field = "lvl1")
fish.net2 <- partic.coef(fish.net2, 
                         db = fishdb1, site.field = "Basin",
                         species.field = "Species",
                         cluster.field = "lvl1")
rbPal <- colorRampPalette(c('blue','red'))
fish.net1$pcCol <- rbPal(10)[as.numeric(cut(fish.net1$participation.coef, breaks = 10))]
basins@data$pcCol1 <- fish.net1$pcCol[match(basins@data$BASIN, fish.net1$Name)]
fish.net2$pcCol <- rbPal(10)[as.numeric(cut(fish.net2$participation.coef, breaks = 10))]
basins@data$pcCol2 <- fish.net2$pcCol[match(basins@data$BASIN, fish.net2$Name)]

op <- par(mfrow = c(2, 1), mar = c(1.1, 1.1, 1.1, 1.1))
plot(basins, col = basins@data$pcCol1)
plot(basins, col = basins@data$pcCol2)
par(op)


