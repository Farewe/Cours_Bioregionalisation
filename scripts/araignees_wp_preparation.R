#### Préparation des données araignées de l'Ouest-Paléarctique pour le TP ####
# Boris Leroy - UMR BOREA, MNHN
#
# 1. Téléchargement des occurrences d'araignées (Araneae) depuis GBIF
# 2. Nettoyage (rang espèce, années, incertitude des coordonnées)
# 3. Maillage sur une grille équivalente en surface (LAEA Europe, EPSG:3035)
# 4. Export au format attendu par le package bioregion :
#    - une matrice sites x espèces (comme bioregion::vegemat)
#    - un objet sf des cellules de la grille (comme bioregion::vegesf)
#
# Prérequis : un compte GBIF (https://www.gbif.org), avec les identifiants
# renseignés dans ~/.Renviron (puis redémarrer R) :
#   GBIF_USER=mon_identifiant
#   GBIF_PWD=mon_mot_de_passe
#   GBIF_EMAIL=mon@email.fr
# usethis::edit_r_environ() ouvre directement le fichier.

library(rgbif)
library(sf)
library(terra)
library(data.table)
library(rnaturalearth)
library(bioregion)

#### Paramètres ####

# Résolution de la grille en km (25, 50, 100...)
resolution_km <- 100

# Filtres d'occurrences
annee_min <- 1970
incertitude_max_m <- 10000 # les incertitudes non renseignées (NA) sont conservées

# Emprise de l'Ouest-Paléarctique (lon lat, sens anti-horaire, exigé par GBIF) :
# Macaronésie et Maghreb jusqu'à ~25°N, Moyen-Orient jusqu'au Zagros,
# Caspienne, fleuve Oural et monts Oural à l'est, Arctique (Svalbard,
# Nouvelle-Zemble) au nord. Le Groenland (néarctique) est exclu plus bas.
wp_wkt <- "POLYGON((-32 25, 48 25, 48 38, 53 47, 60 52, 60 68, 70 82, -32 82, -32 25))"

crs_grille <- "EPSG:3035" # Lambert Azimuthal Equal Area Europe

dossier_gbif <- "./data/gbif_raw"   # téléchargement brut (volumineux, non versionné)
dossier_sortie <- "./data/araignees"

dir.create(dossier_gbif, showWarnings = FALSE, recursive = TRUE)
dir.create(dossier_sortie, showWarnings = FALSE, recursive = TRUE)


#### 1. Téléchargement GBIF ####

# Clé taxonomique de l'ordre Araneae dans le référentiel GBIF (1496)
araneae_key <- name_backbone("Araneae", rank = "order")$usageKey

# La clé du téléchargement est sauvegardée pour ne pas relancer la requête
# (plusieurs minutes à quelques heures côté GBIF) à chaque exécution.
fichier_cle <- file.path(dossier_gbif, "download_key.txt")

if (file.exists(fichier_cle)) {
  dl_key <- readLines(fichier_cle)
} else {
  dl_key <- occ_download(
    pred("taxonKey", araneae_key),
    pred_within(wp_wkt),
    pred("hasCoordinate", TRUE),
    pred("hasGeospatialIssue", FALSE),
    pred("occurrenceStatus", "PRESENT"),
    pred_gte("year", annee_min),
    pred_or(pred_lte("coordinateUncertaintyInMeters", incertitude_max_m),
            pred_isnull("coordinateUncertaintyInMeters")),
    pred_not(pred_in("basisOfRecord",
                     c("FOSSIL_SPECIMEN", "LIVING_SPECIMEN"))),
    format = "SIMPLE_CSV"
  )
  writeLines(dl_key, fichier_cle)
}

occ_download_wait(dl_key)
fichier_zip <- occ_download_get(dl_key, path = dossier_gbif)

# Citation du jeu de données (DOI) à reporter dans le TP
print(gbif_citation(occ_download_meta(dl_key))$download)


#### 2. Lecture et nettoyage ####

colonnes <- c("species", "taxonRank", "decimalLongitude", "decimalLatitude",
              "coordinateUncertaintyInMeters", "year", "basisOfRecord",
              "countryCode")

occ <- occ_download_import(fichier_zip, select = colonnes)
setDT(occ)
nrow(occ)

# Rang espèce uniquement : les sous-espèces sont regroupées dans leur espèce
# (colonne species), les identifications au genre ou à la famille sont retirées
occ <- occ[taxonRank %in% c("SPECIES", "SUBSPECIES", "VARIETY", "FORM") &
             !is.na(species) & species != ""]

# Groenland exclu (néarctique)
occ <- occ[is.na(countryCode) | countryCode != "GL"]

nrow(occ)
uniqueN(occ$species)


#### 3. Grille équivalente en surface ####

# Masque terrestre de l'Ouest-Paléarctique
sf_use_s2(FALSE)
wp <- st_as_sfc(wp_wkt, crs = "EPSG:4326")
terres <- ne_countries(scale = 50, returnclass = "sf")
terres <- terres[terres$iso_a2_eh != "GL", ]
terres <- st_union(st_make_valid(terres))
wp_terres <- st_transform(st_intersection(terres, wp), crs_grille)
sf_use_s2(TRUE)

# Grille alignée sur l'origine de la projection (comme la grille de
# référence de l'Agence Européenne de l'Environnement)
res_m <- resolution_km * 1000
emprise <- ext(vect(wp_terres))
emprise <- ext(floor(emprise$xmin / res_m) * res_m,
               ceiling(emprise$xmax / res_m) * res_m,
               floor(emprise$ymin / res_m) * res_m,
               ceiling(emprise$ymax / res_m) * res_m)
grille <- rast(emprise, resolution = res_m, crs = crs_grille)

# Cellules terrestres : toute cellule touchée par une terre émergée
cellules_terre <- rasterize(vect(wp_terres), grille, touches = TRUE)
cellules_terre <- which(!is.na(values(cellules_terre)))

# Projection des occurrences et attribution d'une cellule
xy <- sf_project(from = "EPSG:4326", to = crs_grille,
                 pts = as.matrix(occ[, .(decimalLongitude, decimalLatitude)]))
occ[, Site := cellFromXY(grille, xy)]

# Les occurrences en mer ou hors emprise sont retirées
occ <- occ[Site %in% cellules_terre]
nrow(occ)

#### 4. Matrice sites x espèces ####

# Nombre d'occurrences par espèce et par cellule (même principe que vegemat)
reseau <- occ[, .(n_occ = .N), by = .(Site, species)]
setorder(reseau, Site, species)
reseau[, Site := as.character(Site)]

comat <- net_to_mat(as.data.frame(reseau), weight = TRUE)

# Transformation en présence-absence
comat[comat > 0] <- 1
storage.mode(comat) <- "integer"
names(dimnames(comat)) <- c("Site", "Species")

dim(comat)

# Polygones des cellules échantillonnées (même principe que vegesf), avec
# l'effort d'échantillonnage pour discuter du biais dans le TP
sites_ids <- as.integer(rownames(comat))
grille_ids <- grille
values(grille_ids) <- NA
grille_ids[sites_ids] <- sites_ids
sites_sf <- st_as_sf(as.polygons(grille_ids, aggregate = FALSE, na.rm = TRUE))
names(sites_sf)[1] <- "Site"
sites_sf$Site <- as.character(sites_sf$Site)
sites_sf <- sites_sf[match(rownames(comat), sites_sf$Site), ]
sites_sf$richesse <- rowSums(comat)
n_occ_site <- reseau[, .(n_occ = sum(n_occ)), by = Site]
sites_sf$n_occurrences <- n_occ_site$n_occ[match(sites_sf$Site,
                                                 n_occ_site$Site)]
sites_sf <- sites_sf[, c("Site", "richesse", "n_occurrences")]
rownames(sites_sf) <- NULL

summary(sites_sf$richesse)


#### 5. Export ####

suffixe <- paste0("wp_", resolution_km, "km")

saveRDS(comat,
        file.path(dossier_sortie, paste0("araignees_", suffixe, "_comat.RDS")))
saveRDS(sites_sf,
        file.path(dossier_sortie, paste0("araignees_", suffixe, "_sites.RDS")))
writeLines(c(paste("GBIF download key:", dl_key),
             gbif_citation(occ_download_meta(dl_key))$download),
           file.path(dossier_sortie, "citation_gbif.txt"))
