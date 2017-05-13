library(maptools)
library(betapart)
library(recluster)
library(dendextend)


# Shapefiles
basins <- readShapeSpatial("./data/data_cours/basin2013_simplif")

# Occurrence databases
load("./data/data_cours/fishdb1.RData")
load("./data/data_cours/fishdb2.RData")
fishdb1 <- as.matrix(table(fishdb1$Basin, fishdb1$Species))
fishdb2 <- as.matrix(table(fishdb2$Basin, fishdb2$Species))


#### Base 1 ####
# 1.
fishdb1 <- fishdb1[-which(rowSums(fishdb1) <= 20), ] # Retrait sites à faible richesse


# 2.
dist.fish1 <- beta.pair(fishdb1, 
                        index.family ="sorensen")


# 3.
fish.nmds1 <- metaMDS(dist.fish1$beta.sim, center=TRUE)

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



# 4.
tree.fish1 <- recluster.cons(dist.fish1$beta.sim, p = 0.5)
plot(tree.fish1$cons,
     tip.color = rgb(red = col.fish1[match(tree.fish1$cons$tip.label, rownames(col.fish1)), 3],
                     green = col.fish1[match(tree.fish1$cons$tip.label, rownames(col.fish1)), 4],
                     blue = col.fish1[match(tree.fish1$cons$tip.label, rownames(col.fish1)), 5],
                     maxColorValue = 255))

# 5.
coph.fish1 <- cophenetic(as.hclust(tree.fish1$cons))
coph.fish1 <- as.matrix(coph.fish1)
coph.fish1 <- coph.fish1[match(attr(dist.fish1$beta.sim, "Labels"), 
                               rownames(coph.fish1)),
                         match(attr(dist.fish1$beta.sim, "Labels"), 
                               colnames(coph.fish1))]
dist.fish1.matrix <- as.matrix(dist.fish1$beta.sim)

cor(dist.fish1.matrix[lower.tri(dist.fish1.matrix)],
    coph.fish1[lower.tri(coph.fish1)])


# 6.
fish.cut1 <- recluster.expl.diss(tree.fish1$cons, dist.fish1$beta.sim, maxcl = 50)

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

# Attention à remettre les noms dans l'ordre alphabétique :
clusters.fish1 <- clusters.fish1[order(names(clusters.fish1))]


# 7.
groupcol.fish1 <- recluster.group.col(col.fish1[1:fish.nmds1$nobj, ], 
                                      clusters.fish1)$all
groupcol.fish1 <- rbind(groupcol.fish1,
                        matrix(data = c(NA, NA, 200, 200, 200),
                               nr = length(which(!(basins@data$BASIN %in% rownames(groupcol.fish1)))),
                               nc = 5, byrow = T,
                               dimnames = list(basins@data$BASIN[which(!(basins@data$BASIN %in% rownames(groupcol.fish1)))])))

# 8.
basins@data$cluster <- NA
basins@data$cluster[match(names(clusters.fish1), basins@data$BASIN)] <- clusters.fish1
basins@data[, c("nmdsx1", "nmdsy1", "r1", "g1", "b1")] <-  groupcol.fish1[match(basins@data$BASIN, rownames(groupcol.fish1)), ]

plot(basins, col = rgb(red = basins@data$r1,
                       green = basins@data$g1,
                       blue = basins@data$b1,
                       maxColorValue = 255))


#### Base 2 ####
# 9.
fishdb2 <- fishdb2[rownames(fishdb2) %in% rownames(fishdb1), ] 

# Distance matrix
dist.fish2 <- beta.pair(fishdb2, 
                        index.family ="sorensen")


# NMDS
fish.nmds2 <- metaMDS(dist.fish2$beta.sim, center=TRUE)

col.fish2 <- recluster.col(fish.nmds2$points)

col.fish2 <- rbind(col.fish2,
                   matrix(data = c(NA, NA, 200, 200, 200),
                          nr = length(which(!(basins@data$BASIN %in% rownames(col.fish2)))),
                          nc = 5, byrow = T,
                          dimnames = list(basins@data$BASIN[which(!(basins@data$BASIN %in% rownames(col.fish2)))])))

op <- par(mfrow = c(2, 1), mar = c(4.1, 4.1, 1.1, 1.1))
recluster.plot.col(col.fish2[1:fish.nmds2$nobj, ])
plot(basins, col = rgb(red = col.fish2[match(basins@data$BASIN, rownames(col.fish2)), 3],
                       green = col.fish2[match(basins@data$BASIN, rownames(col.fish2)), 4],
                       blue = col.fish2[match(basins@data$BASIN, rownames(col.fish2)), 5],
                       maxColorValue = 255))

par(op)

# Hierarchical clustering
# tree.fish1 <- hclust(dist.fish1$beta.sim, method = "average")
# plot(tree.fish1)

# Consensus tree
tree.fish2 <- recluster.cons(dist.fish2$beta.sim, p = 0.5)
# boot.fish1 <- recluster.boot(tree.fish1$cons, dist.fish1$beta.sim, tr = 10, boot = 50)

op <- par(mfrow = c(2, 1), mar = c(1.1, 1.1, 1.1, 1.1))
plot(tree.fish1$cons,
     tip.color = rgb(red = col.fish1[match(tree.fish1$cons$tip.label, rownames(col.fish1)), 3],
                     green = col.fish1[match(tree.fish1$cons$tip.label, rownames(col.fish1)), 4],
                     blue = col.fish1[match(tree.fish1$cons$tip.label, rownames(col.fish1)), 5],
                     maxColorValue = 255))

plot(tree.fish2$cons,
     tip.color = rgb(red = col.fish2[match(tree.fish2$cons$tip.label, rownames(col.fish2)), 3],
                     green = col.fish2[match(tree.fish2$cons$tip.label, rownames(col.fish2)), 4],
                     blue = col.fish2[match(tree.fish2$cons$tip.label, rownames(col.fish2)), 5],
                     maxColorValue = 255))
par(op)

# Cut tree at the same height as before
clusters.fish2 <- cutree(tree.fish2$cons, h = h)
# Attention à remettre les noms dans l'ordre alphabétique :
clusters.fish2 <- clusters.fish2[order(names(clusters.fish2))]

groupcol.fish2 <- recluster.group.col(col.fish2[1:fish.nmds2$nobj, ], 
                                      clusters.fish2)$all

groupcol.fish2 <- rbind(groupcol.fish2,
                        matrix(data = c(NA, NA, 200, 200, 200),
                               nr = length(which(!(basins@data$BASIN %in% rownames(groupcol.fish2)))),
                               nc = 5, byrow = T,
                               dimnames = list(basins@data$BASIN[which(!(basins@data$BASIN %in% rownames(groupcol.fish2)))])))

basins@data$cluster <- NA
basins@data$cluster[match(names(clusters.fish2), basins@data$BASIN)] <- clusters.fish2
basins@data[, c("nmdsx2", "nmdsy2", "r2", "g2", "b2")] <-  groupcol.fish2[match(basins@data$BASIN, rownames(groupcol.fish2)), ]



op <- par(mfrow = c(2, 1), mar = c(1.1, 1.1, 1.1, 1.1))
plot(basins, col = rgb(red = basins@data$r1,
                       green = basins@data$g1,
                       blue = basins@data$b1,
                       maxColorValue = 255))

plot(basins, col = rgb(red = basins@data$r2,
                       green = basins@data$g2,
                       blue = basins@data$b2,
                       maxColorValue = 255))

par(op)





