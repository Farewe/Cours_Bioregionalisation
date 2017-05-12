library(maptools)
library(betapart)
library(recluster)
library(dendextend)

# Occurrence databases
load("./data/data_cours/fishdb1.RData")
contin.fish1 <- as.matrix(table(fishdb1$Basin, fishdb1$Species))
contin.fish1 <- contin.fish1[-which(rowSums(contin.fish1) == 0), ]
contin.fish1 <- contin.fish1[-which(rowSums(contin.fish1) <= 10), ] # Retrait sites à faible richesse
# contin.fish1 <- contin.fish1[-which(rownames(contin.fish1) == "Hiva.Oa3"), ] # Retrait site éloigné 



# Shapefiles
basins <- readShapeSpatial("./data/data_cours/basin2013_simplif")

# Distance matrix
dist.fish1 <- beta.pair(contin.fish1, 
                        index.family ="sorensen")


# NMDS
fish.nmds1 <- metaMDS(dist.fish1$beta.sim, center=TRUE, try = 5, max.try = 5)

col.fish1 <- recluster.col(fish.nmds1$points)

col.fish1 <- rbind(col.fish1,
                  matrix(data = c(NA, NA, 200, 200, 200),
                         nr = length(which(!(basins@data$BASIN %in% rownames(col.fish1)))),
                         nc = 5, byrow = T,
                         dimnames = list(basins@data$BASIN[which(!(basins@data$BASIN %in% rownames(col.fish1)))])))

op <- par(mfrow = c(2, 1), mar = c(4.1, 4.1, 1.1, 1.1))
recluster.plot.col(col.fish1[1:fish.nmds1$nobj, ])
plot(basins, col = rgb(red = col.fish1[match(basins@data$BASIN, rownames(col.fish1)), 3],
                       green = col.fish1[match(basins@data$BASIN, rownames(col.fish1)), 4],
                       blue = col.fish1[match(basins@data$BASIN, rownames(col.fish1)), 5],
                       maxColorValue = 255))

par(op)

# Hierarchical clustering
# tree.fish1 <- hclust(dist.fish1$beta.sim, method = "average")
# plot(tree.fish1)

# Consensus tree
tree.fish1 <- recluster.cons(dist.fish1$beta.sim, p = 0.5)
# boot.fish1 <- recluster.boot(tree.fish1$cons, dist.fish1$beta.sim, tr = 10, boot = 50)
plot(tree.fish1$cons,
     tip.color = rgb(red = col.fish1[match(tree.fish1$cons$tip.label, rownames(col.fish1)), 3],
                     green = col.fish1[match(tree.fish1$cons$tip.label, rownames(col.fish1)), 4],
                     blue = col.fish1[match(tree.fish1$cons$tip.label, rownames(col.fish1)), 5],
                     maxColorValue = 255))

# Corrélation entre distance cophénétique et distance de base
coph.fish1 <- cophenetic(tree.fish1$cons)
cor(dist.fish1$beta.sim, coph.fish1)



# Test different cuts
fish.cut1 <- recluster.expl.diss(tree.fish1$cons, dist.fish1$beta.sim, maxcl = 50)
# library(phytools)
# fish.cut1 <- explore.div(tree.fish1$cons, dist.fish1$beta.sim, maxcl = 40)

plot(fish.cut1$expl.div ~ fish.cut1$nclust, type = "l", lwd = 2, ylim = c(0, 1),
     bty = "l", xlab = "Nombre de clusters", ylab = "% de dissimilarité", las = 1)
abline(h = 0.9, lty = 3)
abline(h = 0.95, lty = 3)
abline(h = 0.99, lty = 3)
abline(h = 0.999, lty = 3)
nclust1 <- fish.cut1$nclust[which(fish.cut1$expl.div > 0.9)][1]
abline(v = nclust1, lty = 2)

k <- 0
h <- 1
while(k < nclust1)
{
  h <- h - .01
  clusters.fish1 <- cutree(tree.fish1$cons, h = h)
  k <- max(clusters.fish1)
}


# Attribute colors to sites according to the selected cutoff
groupcol.fish1 <- recluster.group.col(col.fish1[1:fish.nmds1$nobj, ], 
                                      fish.cut1$matrix[, which(fish.cut1$expl.div > 0.95)[1]])$all
groupcol.fish1 <- rbind(groupcol.fish1,
                        matrix(data = c(NA, NA, 200, 200, 200),
                               nr = length(which(!(basins@data$BASIN %in% rownames(groupcol.fish1)))),
                               nc = 5, byrow = T,
                               dimnames = list(basins@data$BASIN[which(!(basins@data$BASIN %in% rownames(groupcol.fish1)))])))

basins@data$cluster <- NA
basins@data$cluster[match(names(clusters.fish1), basins@data$BASIN)] <- clusters.fish1
basins@data[, c("nmdsx", "nmdsy", "r", "g", "b")] <-  groupcol.fish1[match(basins@data$BASIN, rownames(groupcol.fish1)), ]

plot(basins, col = rgb(red = basins@data$r,
                       green = basins@data$g,
                       blue = basins@data$b,
                       maxColorValue = 255))


