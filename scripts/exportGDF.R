exportGDF <- function (network, clusters, filename, weight = FALSE,
                       hex2rgb = TRUE,
                       color.field = NULL,
                       directed = FALSE)
  # , network = NULL, site.field = colnames(db)[1], species.field = colnames(db)[2], 
  #         filename, color.field = NULL, abundance.field = NULL, hex2rgb = TRUE, 
  #         directed = FALSE, additional.fields = NULL) 
{
  scipen <- options()$scipen
  options(scipen = 999999)
  
  if (any(is.na(network))) 
    network[is.na(network)] <- 0
  if (length(color.field)  & hex2rgb) {
    clusters[, color.field] <- 
      apply(grDevices::col2rgb(tolower(clusters[, color.field])),
            2, paste, collapse = ",")
    clusters[, color.field] <- 
      paste0("'",     clusters[, color.field], 
             "'")
  }
  species.table <- clusters[which(attributes(clusters)$node_type == "species"), ] 
  sites.table <- clusters[which(attributes(clusters)$node_type == "site"), ]  
  
  
  species.table$id <- 1:nrow(species.table)
  sites.table$id <- (nrow(species.table) + 1):(nrow(sites.table) + nrow(species.table))
  
  links <- data.frame(
    from = species.table$id[match(network[, "Species"],
                                  species.table$ID)], 
    to = sites.table$id[match(network[, "Basin"], 
                              sites.table$ID)], 
    weight = ifelse(rep(weight, nrow(network)),
                    network[, 3], rep(1, nrow(network))))
  
  links[, c(1:2)] <- sapply(links[, c(1:2)], as.integer)
  if (directed) {
    links$direction <- rep("true", nrow(links))
  }
  cat(paste("nodedef>name VARCHAR,label VARCHAR,", 
            paste(colnames(clusters)[grep("K_", colnames(clusters))], " VARCHAR", 
                  sep = "", collapse = ","),
            ifelse(length(color.field), 
                   ifelse(hex2rgb, ",color VARCHAR", ",ccolor VARCHAR"), 
                   ""),
            "\n",
            paste0(paste(species.table$id,
                         species.table$ID,
                         apply(species.table[, grep("K_", colnames(species.table))], 
                               1, paste, collapse = ","), sep = ","), 
                   ifelse(rep(length(color.field), 
                              nrow(species.table)), 
                          paste0(",", 
                                 species.table[, color.field]), 
                          ""),
                   collapse = "\n"), "\n", 
            paste0(paste(sites.table$id, 
                         sites.table$Name, 
                         apply(sites.table[, grep("K_", 
                                                  colnames(sites.table))], 
                               1, paste, collapse = ","), 
                         sep = ","), 
                   ifelse(rep(length(color.field), nrow(sites.table)), 
                          paste0(",", 
                                 sites.table[, color.field]),
                          ""), collapse = "\n"), 
            "\n", paste0("edgedef>node1 VARCHAR,node2 VARCHAR,weight INTEGER", 
                         ifelse(directed, ",directed BOOLEAN", "")), "\n", 
            paste(apply(links, 1, paste, collapse = ","), collapse = "\n"), 
            sep = ""), file = filename)
  
  # 
  # if (!is.null(network)) {
  #   network[, grep("lvl", colnames(network))] <- apply(network[, 
  #                                                              grep("lvl", colnames(network))], 2, as.character)
  #   if (any(is.na(network))) 
  #     network[is.na(network)] <- 0
  #   if (length(color.field) & hex2rgb) {
  #     network[, color.field] <- apply(col2rgb(tolower(network[, 
  #                                                             color.field])), 2, paste, collapse = ",")
  #     network[, color.field] <- paste0("'", network[, color.field], 
  #                                      "'")
  #   }
  #   species.table <- getSpeciesTable(db = db, network = network, 
  #                                    species.field = species.field)
  #   sites.table <- getSiteTable(db = db, network = network, 
  #                               site.field = site.field)
  #   links <- data.frame(from = species.table$id[match(db[, 
  #                                                        species.field], species.table$Name)], to = sites.table$id[match(db[, 
  #                                                                                                                           site.field], sites.table$Name)], weight = ifelse(rep(length(abundance.field), 
  #                                                                                                                                                                                nrow(db)), db[, abundance.field], rep(1, nrow(db))))
  #   links[, c(1:2)] <- sapply(links[, c(1:2)], as.integer)
  #   if (directed) {
  #     links$direction <- rep("true", nrow(links))
  #   }
  #   cat(paste("nodedef>name VARCHAR,label VARCHAR,", 
  #             paste("lvl", 
  #                   1:length(grep("lvl", colnames(network))), " VARCHAR", 
  #                   sep = "", collapse = ","),
  #             ifelse(length(color.field), 
  #                    ifelse(hex2rgb, ",color VARCHAR", ",ccolor VARCHAR"), 
  #                    ""),
  #             "\n",
  #             paste0(paste(species.table$id,
  #                          species.table$Name, 
  #                          
  #                          apply(species.table[, grep("lvl", colnames(species.table))], 
  #                                1, paste, collapse = ","), sep = ","), 
  #                    ifelse(rep(length(color.field), 
  #                               nrow(species.table)), paste0(",", 
  #                                                            species.table[, 
  #                                                                          color.field]), ""),
  #                    collapse = "\n"), "\n", 
  #             paste0(paste(sites.table$id, 
  #                          sites.table$Name, apply(sites.table[, grep("lvl", 
  #                                                                     colnames(sites.table))], 1, paste, collapse = ","), 
  #                          sep = ","), ifelse(rep(length(color.field), nrow(sites.table)), 
  #                                             paste0(",", sites.table[, color.field]), ""), collapse = "\n"), 
  #             "\n", paste0("edgedef>node1 VARCHAR,node2 VARCHAR,weight INTEGER", 
  #                          ifelse(directed, ",directed BOOLEAN", "")), "\n", 
  #             paste(apply(links, 1, paste, collapse = ","), collapse = "\n"), 
  #             sep = ""), file = filename)
  # }
  # else {
  #   if (directed) {
  #     db$direction <- rep("true", nrow(db))
  #   }
  #   else {
  #     direction <- NULL
  #   }
  #   if (length(color.field)) {
  #     colors <- cbind(unique(db[, c(site.field, color.field)]), 
  #                     deparse.level = 0)
  #     colnames(colors) <- c("node", "color")
  #     if (any(!(db[, species.field] %in% colors[, 1]))) {
  #       colors <- rbind(colors, cbind(node = unique(as.character(db[which(!(db[, 
  #                                                                              species.field] %in% colors[, 1])), species.field])), 
  #                                     color = grey(0.5)))
  #     }
  #     if (hex2rgb) {
  #       colors[, "color"] <- apply(col2rgb(tolower(colors[, 
  #                                                         "color"])), 2, paste, collapse = ",")
  #       colors[, "color"] <- paste0("'", colors[, "color"], 
  #                                   "'")
  #     }
  #   }
  #   if (!length(abundance.field)) {
  #     abundance.field <- "weight"
  #     db$weight <- rep(1, nrow(db))
  #   }
  #   cat(paste(paste0("nodedef>name VARCHAR", ifelse(length(color.field), 
  #                                                   ",color VARCHAR\n", "\n")), paste(ifelse(rep(length(color.field), 
  #                                                                                                length(unique(c(levels(as.factor(db[, site.field])), 
  #                                                                                                                levels(as.factor(db[, species.field])))))), apply(colors, 
  #                                                                                                                                                                  1, paste, collapse = ","), unique(c(levels(as.factor(db[, 
  #                                                                                                                                                                                                                          site.field])), levels(as.factor(db[, species.field]))))), 
  #                                                                                     collapse = "\n"), "\n", paste0("edgedef>node1 VARHCHAR,node2 VARCHAR,weight INTEGER", 
  #                                                                                                                    ifelse(directed, ",directed BOOLEAN", ""), "\n"), 
  #             paste(apply(db[, c(site.field, species.field, abundance.field, 
  #                                direction)], 1, paste, collapse = ","), collapse = "\n"), 
  #             sep = ""), file = filename)
  # }
  options(scipen = scipen)
}