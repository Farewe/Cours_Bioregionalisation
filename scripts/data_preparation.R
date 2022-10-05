#### species data ####

a <- read.csv('./data/database/Occurences.csv', sep = ";")

base2016 <- read.table("j:/r/projects/FishBiogeography/data/BasePoisson102016_2.csv", h = T, sep =  ";")

base2019 <- read.table("d:/r/Projects/AnthropoceneFish/data/Occurrence_Table_29_08_2019.csv", h = T, sep =  ";")

fish.all <- read.csv("d:/r/Projects/AnthropoceneFish/data/Occurrence_Table_29_08_2019.csv", h = T, sep =  ";")

fish.info <- unique(base2016[, c("FB_REF2015", "fresh", "brack", "salt", "AnaCat")])


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

fish.all <- readRDS("d:/r/Projects/AnthropoceneFish/data/fish_anthropocene.rds")

cont <- table(fish.all$X1.Basin.Name,
              fish.all$X6.Fishbase.Valid.Species.Name)

richness <- rowSums(cont)

low.rich.basins <- names(richness)[which(richness <= 20)]

fish.all <- as.data.frame(fish.all[-which(fish.all$X1.Basin.Name %in% low.rich.basins), ])

fish.all <- fish.all[, c("X1.Basin.Name", "X6.Fishbase.Valid.Species.Name", "X3.Native.Exotic.Status")]
fish.all$X6.Fishbase.Valid.Species.Name <- as.factor(fish.all$X6.Fishbase.Valid.Species.Name)
levels(fish.all$X6.Fishbase.Valid.Species.Name) <- paste0("sp", 1:nlevels(fish.all$X6.Fishbase.Valid.Species.Name))
fish.all$X1.Basin.Name <- as.factor(fish.all$X1.Basin.Name)

names(fishdb1)
fishdb1 <- droplevels(fish.all[which(fish.all$X3.Native.Exotic.Status == "native"), c(1, 2)])
fishdb2 <- droplevels(fish.all[, c(1, 2)])

names(fishdb1) <- names(fishdb2) <- c("Basin", "Species")

save(fishdb1, file = "./data/data_cours/fishdb1.RData")
save(fishdb2, file = "./data/data_cours/fishdb2.RData")





#### shapefile data ####
library(rgdal)
basinshp <- readOGR("d:/r/Projects/AnthropoceneFish/data/sig data/basinshpwithbaikal")

levels(fishdb1$Basin)[!(levels(fishdb1$Basin) %in% unique(basinshp@data$BasinName))]
basinshp <- basinshp[-which(!(basinshp$BasinName %in% levels(fishdb1$Basin))), ]
unique(basinshp@data$BasinName)[!(unique(basinshp@data$BasinName) %in% levels(fishdb1$Basin))]

writeOGR(basinshp, "./data/data_cours", "basin_simplif", driver="ESRI Shapefile")

library(sf)
basinshp_simp <- st_make_valid(st_as_sf(basinshp))

basinshp_simp <- st_simplify(basinshp_simp,
                             dTolerance = 20000)
plot(basinshp_simp[1])

st_write(basinshp_simp, "./data/data_cours",
         "basin_simplif", driver = "ESRI Shapefile")

writeOGR(basinshp_simp, "./data/data_cours", "basin_simplif", driver="ESRI Shapefile")

