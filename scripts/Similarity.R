
library(betapart)
library(recluster)
library(dendextend)
library(vegan)
library(rnaturalearth)
devtools::source_url("https://raw.githubusercontent.com/Farewe/Cours_Bioregionalisation/master/scripts/modified_recluster_tree_function.R")

wm <- getMap(resolution = "low")

# Occurrence database
subtidaldb <- readRDS("./data/invertebres_benthiques.RDS")


# Site position
sites <- readRDS("./data/sites_invertebresbenthiques.RDS")

#### Base 1 ####
# 1.
# fishdb1 <- fishdb1[-which(rowSums(fishdb1) <= 20), ] # Retrait sites à faible richesse


# 2.
dist.subtidal <- beta.pair(t(subtidaldb), 
                        index.family ="sorensen")

# 3.
subtidal.nmds <- metaMDS(dist.subtidal$beta.sim, center=TRUE)

col.subtidal <- recluster.col(subtidal.nmds$points)

colnames(col.subtidal) <- c("nmdsx", "nmdsy", "nmdsred", "nmdsgreen", "nmdsblue")

sites <- data.frame(sites, 
                    col.subtidal[match(sites$station_id,
                                       rownames(col.subtidal)), ])


op <- par(mfrow = c(2, 1), mar = c(4.1, 4.1, 1.1, 1.1))
recluster.plot.col(col.subtidal[1:subtidal.nmds$nobj, ])

plot(sites$Y ~ sites$X,
     col = rgb(red = sites$nmdsred,
               green = sites$nmdsgreen,
               blue = sites$nmdsblue,
               maxColorValue = 255))
wm <- ne_coastline(scale = 50, returnclass = "sp")

plot(wm, add = TRUE)

par(op)



# 4.
tree.subtidal <- modified.recluster.cons(dist.subtidal$beta.sim, tr = 100)




# 5.
coph.subtidal <- as.matrix(cophenetic(as.hclust(tree.subtidal$cons)))
# Pour corriger l'ordre des sites
coph.subtidal <- coph.subtidal[match(attr(dist.subtidal$beta.sim, "Labels"), 
                                     rownames(coph.subtidal)),
                               match(attr(dist.subtidal$beta.sim, "Labels"), 
                                     colnames(coph.subtidal))]
dist.subtidal.matrix <- as.matrix(dist.subtidal$beta.sim)

cor(dist.subtidal.matrix[lower.tri(dist.subtidal.matrix)],
    coph.subtidal[lower.tri(coph.subtidal)])

max(tree.subtidal$cophcor)


chosen.tree <- tree.subtidal$trees[[which(tree.subtidal$cophcor == max(tree.subtidal$cophcor))]]



plot(chosen.tree,
     tip.color = rgb(red = col.subtidal[match(chosen.tree$tip.label, rownames(col.subtidal)), 3],
                     green = col.subtidal[match(chosen.tree$tip.label, rownames(col.subtidal)), 4],
                     blue = col.subtidal[match(chosen.tree$tip.label, rownames(col.subtidal)), 5],
                     maxColorValue = 255))

axis(1, 
     at = seq(max(cophenetic(as.hclust(tree.subtidal$cons))), 0, 
              by = -round(max(cophenetic(as.hclust(tree.subtidal$cons)))/4, 2)) / 2,
     labels = seq(0, max(cophenetic(as.hclust(tree.subtidal$cons))), 
                  by = round(max(cophenetic(as.hclust(tree.subtidal$cons)))/4, 2)))



# 6.
# Méthode de Holt et al. 2013
cut.subtidal <- recluster.expl.diss(chosen.tree, dist.subtidal$beta.sim, maxcl = 50)

plot(cut.subtidal$expl.div ~ cut.subtidal$nclust, type = "l", lwd = 2, ylim = c(0, 1),
     bty = "l", xlab = "Nombre de clusters", ylab = "% de dissimilarité", las = 1)
abline(h = 0.9, lty = 3)
abline(h = 0.95, lty = 3)
abline(h = 0.99, lty = 3)
abline(h = 0.999, lty = 3)
nclust <- cut.subtidal$nclust[which(cut.subtidal$expl.div > 0.9)][1]
abline(v = nclust, lty = 2)


k <- 0
h <- 1
while(k < nclust)
{
  h <- h - .01
  clusters.subtidal <- cutree(chosen.tree, h = h)
  k <- max(clusters.subtidal)
}

h



plot(chosen.tree,
     tip.color = rgb(red = col.subtidal[match(chosen.tree$tip.label, rownames(col.subtidal)), 3],
                     green = col.subtidal[match(chosen.tree$tip.label, rownames(col.subtidal)), 4],
                     blue = col.subtidal[match(chosen.tree$tip.label, rownames(col.subtidal)), 5],
                     maxColorValue = 255))

axis(1, 
     at = seq(max(cophenetic(as.hclust(tree.subtidal$cons))), 0, 
              by = -round(max(cophenetic(as.hclust(tree.subtidal$cons)))/4, 2)) / 2,
     labels = seq(0, max(cophenetic(as.hclust(tree.subtidal$cons))), 
                  by = round(max(cophenetic(as.hclust(tree.subtidal$cons)))/4, 2)))
abline(v = (max(cophenetic(as.hclust(tree.subtidal$cons))) - h) / 2)



# Attention à remettre les noms dans l’ordre alphabétique pour l’étape 7 (car la NMDS a les noms dans l’ordre alphabétique)
clusters.subtidal <- clusters.subtidal[order(names(clusters.subtidal))]


# 7.
groupcol.subtidal <- recluster.group.col(col.subtidal[1:subtidal.nmds$nobj, ], 
                                      clusters.subtidal)$all

# 8.
sites$cluster[match(names(clusters.subtidal), sites$station_id)] <- clusters.subtidal
sites[, c("nmdsx1", "nmdsy1", "r1", "g1", "b1")] <-  
  groupcol.subtidal[match(sites$station_id, rownames(groupcol.subtidal)), ]


plot(sites$Y ~ sites$X, col = rgb(red = sites$r1,
                       green = sites$g1,
                       blue = sites$b1,
                       maxColorValue = 255))
plot(wm, add = TRUE)
