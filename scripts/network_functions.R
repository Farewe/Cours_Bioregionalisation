writePajek <- function(db, site.field = 1, species.field = 2, filename, abundance.field = NULL)
{
  species <- data.frame(sp = levels(db[, species.field]),
                        id = 1:length(levels(db[, species.field])))
  sites <- data.frame(site = levels(db[, site.field]),
                      id = length(levels(db[, species.field])) + 1:length(levels(db[, site.field])))
  
  links <- data.frame(from = species$id[match(db[, species.field], species$sp)], 
                      to = sites$id[match(db[, site.field], sites$site)], 
                      weight = ifelse(rep(length(abundance.field), nrow(db)),
                                      db[, abundance.field],
                                      rep(1, nrow(db))))
                        
  cat(paste("*Vertices ", max(sites$id), "\n",
            paste(species$id, ' "', species$sp, '"', sep = "", collapse = "\n"), "\n",
            paste(sites$id, ' "', sites$site, '"', sep = "", collapse = "\n"), "\n",
            "*Edges\n",
            paste(apply(links, 1, paste, collapse = " "), collapse = "\n"), 
            sep = ""), file = filename)
}

read.infomap.tree <- function(file, max.length = -1, reorder.levels = TRUE)
{
  tree.infomap <- read.table(file, skip = 1, sep = " ", nrows = max.length)
  colnames(tree.infomap) <- c("Groups", "Codelength", "Name", "id")
  tmp1 <- strsplit(as.character(tree.infomap$Groups), ":")
  indx <- sapply(tmp1, length)
  tmp2 <-  as.data.frame(do.call(rbind,lapply(tmp1, `length<-`, max(indx))))
  colnames(tmp2) <- paste("lvl", 1:ncol(tmp2), sep = "") 
  tree.infomap <- data.frame(tree.infomap,
                             tmp2)
  if(reorder.levels)
  {
    tree.infomap[, grep("lvl", colnames(tree.infomap))] <- lapply(tree.infomap[, grep("lvl", colnames(tree.infomap))],
                                                                  function(xx) factor(xx, levels = 1:max(as.numeric(xx), na.rm = T)))
  }
  return(tree.infomap)
}

groupNodesPerCluster <- function(db, site.field = 1, species.field = 2, network, color = FALSE, abundance.field = NULL)
{
  library(plyr)
  db[, c(site.field, species.field)] <- apply(db[, c(site.field, species.field)], 2, as.character)
  db[, site.field] <- network$lvl1[match(db[, site.field], network$Name)]
  db[, species.field] <- network$lvl1[match(db[, species.field], network$Name)]
  db[, c(site.field, species.field)] <- lapply(db[, c(site.field, species.field)], factor)
  if(!is.null(abundance.field))
  {
    cluster.db <- ddply(db, c(site.field, species.field), .fun =
                          function(x){
                            sum(x[, abundance.field])
                                                 }
                        , .progress = "text")
  } else
  {
    cluster.db <- ddply(db, c(site.field, species.field), .fun =
                          function(x){
                            sum(x[, abundance.field])
                          }
                        , .progress = "text")
  }

  if(color)
  {
    cluster.db$color <- apply(col2rgb(tolower(network$color[match(cluster.db[, site.field], network$lvl1)])), 2, paste, collapse = ",")
    cluster.db$color <- paste0("'", cluster.db$color, "'")
  }
  return(cluster.db)
}

attributeColors <- function(network, nb.max.colors = 12, palette = "Paired", lvl = "lvl1")
{
  require(RColorBrewer)
  network$color <-  c(brewer.pal(min(c(max(as.numeric(network[, lvl])), nb.max.colors)), palette), 
                      rep(grey(.5), ifelse(max(as.numeric(network[, lvl])) > nb.max.colors,
                                           max(as.numeric(network[, lvl])) - nb.max.colors, 0)))[as.numeric(network$lvl1)]
  return(network)
}

getSpeciesTable <- function(db, species.field = 2, network)
{
  species <- as.character(network$Name[which(network$Name %in% db[, species.field])])
  species <- network[which(network$Name %in% species), ]
  return(species)
}

getSiteTable <- function(db, site.field = 1, network)
{
  sites <- as.character(network$Name[which(network$Name %in% db[, site.field])])
  sites <- network[which(network$Name %in% sites), ]
  return(sites)
}

create.igraph.from.db <- function(db, network = NULL, site.field = 1, species.field = 2, gp.color = NULL, diff = FALSE)
{
  require(igraph)
  if(!missing(network))
  {
    network[, grep("lvl", colnames(network))] <- apply(network[, grep("lvl", colnames(network))], 2, as.character)
    network[, grep("lvl", colnames(network))] <- apply(network[, grep("lvl", colnames(network))], 2, as.numeric)
    if(any(is.na(network))) network[is.na(network)] <- 0
  }
  species.table <- data.frame(node = levels(db[, species.field]),
                              id = 1:length(levels(db[, species.field])))
  sites.table <- data.frame(node = levels(db[, site.field]),
                            id = length(levels(db[, species.field])) + 1:length(levels(db[, site.field])))
  all.table <- rbind.data.frame(species.table,
                                sites.table)
  links <- data.frame(from = species.table$id[match(db[, species.field], species.table$node)], 
                      to = sites.table$id[match(db[, site.field], sites.table$node)])
  out.graph <- graph.data.frame(links, vertices = all.table[, 2:1], directed = F)
  if(diff)
  {
    V(out.graph)$shape <- c(rep("circle", nrow(species.table)),
                            rep("square", nrow(sites.table)))
  }
  if(length(gp.color))
  {
    V(out.graph)$color <- brewer.pal(max(network[, gp.color]), "Set3")[network[, gp.color][match(all.table$node, 
                                                                                                 network$Name)]]
  }
  return(out.graph)
}

writeGDF <- function(db, network = NULL, site.field = colnames(db)[1], species.field = colnames(db)[2], 
                     filename, color.field = NULL, abundance.field = NULL, hex2rgb = TRUE, directed = FALSE,
                     additional.fields = NULL) # additional.fields to be implemented later
{
  if(!is.null(network)) # When a network is present, then write a gdf file taking into account network levels
  {
    network[, grep("lvl", colnames(network))] <- apply(network[, grep("lvl", colnames(network))], 2, as.character)
    network[, grep("lvl", colnames(network))] <- apply(network[, grep("lvl", colnames(network))], 2, as.numeric)
    if(any(is.na(network))) network[is.na(network)] <- 0
    
    if(length(color.field) & hex2rgb)
    {
      network[, color.field] <- apply(col2rgb(tolower(network[, color.field])), 2, paste, collapse = ",")
      network[, color.field] <- paste0("'", network[, color.field] , "'")
    }
    
    species.table <- getSpeciesTable(db = db, network = network, species.field = species.field)
    sites.table <- getSiteTable(db = db, network = network, site.field = site.field)
    
    links <- data.frame(from = species.table$id[match(db[, species.field], species.table$Name)], # from species
                        to = sites.table$id[match(db[, site.field], sites.table$Name)], # to site
                        weight = ifelse(rep(length(abundance.field), nrow(db)), # if abundance
                                        db[, abundance.field], # paste abundance in weight
                                        rep(1, nrow(db)))) # else set weight = 1
    if(directed)
    {
      links$direction <- rep("true", nrow(links))
    }
    
    cat(paste(# Definition of node fields
      "nodedef>name VARCHAR,label VARCHAR,", 
      paste("lvl", 1:length(grep("lvl", colnames(network))), # indicate level fields
            " INTEGER", sep  = "", 
            collapse = ","), 
      ifelse(length(color.field), # if colors: indicate color field
             ifelse(hex2rgb, # if conversion to gephi format (rgb): name the column 'color', else 'ccolor'
                    ",color VARCHAR",
                    ",ccolor VARCHAR"), 
             ""),
      "\n",
      # Nodes
      ## Species
      paste0(paste(species.table$id, # id
                   species.table$Name, # name
                   apply(species.table[, grep("lvl", colnames(species.table))], 1, paste, collapse = ","), # levels 
                   sep = ","), ifelse(rep(length(color.field), nrow(species.table)), # if colors: paste colors from species table, else do nothing
                                      paste0(",", species.table[, color.field]),
                                      ""), collapse = "\n"),
      "\n", 
      ## Sites
      paste0(paste(sites.table$id, # id
                   sites.table$Name, # name
                   apply(sites.table[, grep("lvl", colnames(sites.table))], 1, paste, collapse = ","), 
                   sep = ","), ifelse(rep(length(color.field), nrow(sites.table)), # if colors: paste colors from site table, else do nothing
                                      paste0(",", sites.table[, color.field]),
                                      ""), collapse = "\n"), 
      "\n",
      # Definition of edge fields
      paste0("edgedef>node1 VARCHAR,node2 VARCHAR,weight INTEGER", 
             ifelse(directed, # if network is directed, indicate it, else do nothing
                    ",directed BOOLEAN",
                    "")),
      "\n",
      # Edge fields: paste everything including weight & direction
      paste(apply(links, 1, paste, collapse = ","), collapse = "\n"), 
      sep = ""), file = filename)
  } else # Else when no network is present, use only the db to create a network and write the gdf file.
  {      # In this case, site.field is the source node field, and species.field is the target node field
    if(directed)
    {
      db$direction <- rep("true", nrow(db))
    } else {direction <- NULL}
    if(length(color.field))
    {
      colors <- cbind(unique(db[, c(site.field, color.field)]), deparse.level = 0)
      colnames(colors) <- c("node", "color")
      if(any(!(db[, species.field] %in% colors[, 1])))
      {
        colors <- rbind(colors,
                        cbind(node = unique(as.character(db[which(!(db[, species.field] %in% colors[, 1])), species.field])),
                              color = grey(.5))) # Attributing grey to nodes with no colour
      }
      if(hex2rgb)
      {
        colors[, "color"] <- apply(col2rgb(tolower(colors[, "color"])), 2, paste, collapse = ",")
        colors[, "color"] <- paste0("'", colors[, "color"], "'")
      }
    }
    if(!length(abundance.field))
    {
      abundance.field <- "weight"
      db$weight <- rep(1, nrow(db))
    }

    
    cat(paste(paste0("nodedef>name VARCHAR", ifelse(length(color.field),
                                                    ",color VARCHAR\n",
                                                    "\n")),
              paste(ifelse(rep(length(color.field), length(unique(c(levels(as.factor(db[, site.field])),
                                                                    levels(as.factor(db[, species.field])))))),
                           apply(colors, 1, paste, collapse = ","),
                           unique(c(levels(as.factor(db[, site.field])),
                                    levels(as.factor(db[, species.field]))))), collapse = "\n"), 
              "\n",
              paste0("edgedef>node1 VARHCHAR,node2 VARCHAR,weight INTEGER",
                     ifelse(directed, ",directed BOOLEAN", ""), "\n"),
              paste(apply(db[, c(site.field, species.field, abundance.field, direction)], 1, paste, collapse = ","), collapse = "\n"), 
              sep = ""), file = filename)
  }
}

writeGDF2 <- function(sharedsptable, network, filename, color.field = NULL,
                      coords)
{
  network[, grep("lvl", colnames(network))] <- apply(network[, grep("lvl", colnames(network))], 2, as.character)
  network[, grep("lvl", colnames(network))] <- apply(network[, grep("lvl", colnames(network))], 2, as.numeric)
  if(any(is.na(network))) network[is.na(network)] <- 0
  
  sharedsptable[upper.tri(sharedsptable, diag = TRUE)] <- 0
  sites.table <- data.frame(site = colnames(sharedsptable),
                            id = 1:ncol(sharedsptable))
  sites.table$x <- coords$x[match(sites.table$site, rownames(coords))]
  sites.table$y <- coords$y[match(sites.table$site, rownames(coords))]

  links <- as.data.frame.table(sharedsptable, stringsAsFactors = FALSE)
  links <- links[-which(links$Freq == 0), ]
  links$Var1 <- sites.table$id[match(links$Var1, sites.table$site)]
  links$Var2 <- sites.table$id[match(links$Var2, sites.table$site)]
 

  if(length(color.field))
  {
    sites.table <- data.frame(sites.table, 
                              network[match(sites.table$site, network$Name), c(which(colnames(network) == "Codelength"),
                                                                               grep("lvl", colnames(network)),
                                                                               which(colnames(network) == color.field))])
    cat(paste("nodedef>name VARCHAR,label VARCHAR,latitude DOUBLE,longitude DOUBLE,", 
              paste("lvl", 1:length(grep("lvl", colnames(network))), " INTEGER", sep  = "", collapse = ","), ",ccolor VARCHAR\n",
              paste(sites.table$id, 
                    sites.table$site, 
                    sites.table$y,
                    sites.table$x,
                    apply(sites.table[, grep("lvl", colnames(sites.table))], 1, paste, collapse = ","), 
                    sep = ",", sites.table[, color.field], collapse = "\n"), "\n",
              "edgedef>node1 VARCHAR,node2 VARCHAR,weight INTEGER\n",
              paste(apply(links, 1, paste, collapse = ","), collapse = "\n"), 
              sep = ""), file = filename)
  } else {
    sites.table <- data.frame(sites.table, 
                              network[match(sites.table$site, network$Name), c(which(colnames(network) == "Codelength"),
                                                                               grep("lvl", colnames(network)))])
    cat(paste("nodedef>name VARCHAR,label VARCHAR,latitude DOUBLE,longitude DOUBLE,", 
              paste("lvl", 1:length(grep("lvl", colnames(network))), " INTEGER", sep  = "", collapse = ","), "\n",
              paste(sites.table$id, 
                    sites.table$site,  
                    sites.table$y,
                    sites.table$x,
                    apply(sites.table[, grep("lvl", colnames(sites.table))], 1, paste, collapse = ","), 
                    sep = ",", collapse = "\n"), "\n",
              "edgedef>node1 VARCHAR,node2 VARCHAR,weight INTEGER\n",
              paste(apply(links, 1, paste, collapse = ","), collapse = "\n"), 
              sep = ""), file = filename)
  }
  
}




remove.low.levels <- function(df, column, limit)
{
  if(length(which(sapply(levels(df[, column]),
                         FUN = function(x, y) length(y[y == x]), 
                         y = df[, column]) < limit)))
  {

    a <- droplevels(df[which(df[, column] %in% 
                               names(sapply(levels(df[, column]), 
                                            FUN = function(x, y) length(y[y == x]), 
                                            y = df[, column]))[which(sapply(levels(df[, column]),
                                                                            FUN = function(x, y) length(y[y == x]), 
                                                                            y = df[, column]) < limit)]), ])
    message(paste(nrow(a), " sites removed, correponding to ", length(levels(a[, column])), " levels.\n",
              "Removed sites: ", paste(levels(a$Name), collapse = ", "), "\n", sep =""))
    df <- df[-which(df[, column] %in% 
                      names(sapply(levels(df[, column]), 
                                   FUN = function(x, y) length(y[y == x]), 
                                   y = df[, column]))[which(sapply(levels(df[, column]),
                                                                     FUN = function(x, y) length(y[y == x]), 
                                                                     y = df[, column]) < limit)]), ]
    df <- droplevels(df)
    return(df)
  } else
  {
    message("No levels below the specified limit\n")
    return(df)
  }
}


plotRadial <- function(network, levels, lvl0 = "world", root.level = NULL,
                       leaf.names = NULL, ...)
{
  require(data.tree)
  require(networkD3)
  
  if(any(!sapply(network[, levels], is.factor)))
  {
    network[, levels] <- lapply(network[, levels], factor)
  }
  network <- droplevels(network)
  
  if(!is.null(root.level))
  {
    if(!(root.level %in% levels(network[, levels[1]])))
    {
      stop("root.level does not exist in your level 1")
    } else
    {
      network <- droplevels(network[which(network[, levels[1]] == root.level), ])
    }
  }
  
  if(!is.null(leaf.names))
  {
    if(!leaf.names %in% colnames(network))
    {
      stop("leaf.names must be the name of one of your dataframe columns")
    } else
    {
      network[, levels] <- lapply(network[, levels], as.character)
      for(i in levels(network[, leaf.names]))
      {
        network[network[, leaf.names] %in% i, levels] <- replace.leaves(i, network, levels, leaf.names)
      }
      network[, levels] <- lapply(network[, levels], factor)
    }
  }
  
  if(length(levels(network[, levels[1]])) == 1)
  {
    network$pathString <- apply(network[, levels], 1, paste, collapse = "/")
  } else
  {
    network$pathString <- paste(lvl0, apply(network[, levels], 1, paste, collapse = "/"), sep = "/")
  }
  
  network$pathString <- gsub("/NA", "", network$pathString)
  
  network.tree <- as.Node(network)
  network.tree.list <- ToListExplicit(network.tree, unname = TRUE)
  
  radialNetwork(network.tree.list, ...)
}

replace.leaves <- function(leaf.name, netwk, lvls, lf.nm)
{
  network.line <- netwk[which(netwk[, lf.nm] == leaf.name), lvls]
  if(any(is.na(network.line)))
  {
    network.line[min(which(is.na(network.line))) - 1] <- leaf.name
  } else
  {
    network.line[length(network.line)] <- leaf.name
  }
  return(network.line)
}

partic.coef <- function(network, db, site.field, species.field,
                        cluster.field = "lvl1")
{
  sites <- as.character(network$Name[which(network$Name %in% db[, site.field])])
  species <- as.character(network$Name[which(network$Name %in% db[, species.field])])
  db$species.cluster <- network[match(db[, species.field], network$Name), cluster.field]
  db$site.cluster <- network[match(db[, site.field], network$Name), cluster.field]
  if(any(sites %in% species) | any(species %in% sites))
  {
    stop("The database does not appear to be a bipartite database")
  }
  network$participation.coef <- NA
  for (site in sites)
  {
    subdb <- db[which(db[, site.field] %in% site), ]
    tmp <- ddply(subdb, "species.cluster", summarise, nlinks = length(species.cluster))
    network$participation.coef[network$Name == site] <- 1 - sum((tmp$nlinks / sum(tmp$nlinks))^2)
  }
  for (sp in species)
  {
    subdb <- db[which(db[, species.field] %in% sp), ]
    tmp <- ddply(subdb, "site.cluster", summarise, nlinks = length(site.cluster))
    network$participation.coef[network$Name == sp] <- 1 - sum((tmp$nlinks / sum(tmp$nlinks))^2)
  }
  return(network)
}

#####################Force Atlas 2 function ########################################################

layout.forceatlas2 <- function(g, iterations = 100, linlog = FALSE, pos = NULL, nohubs = FALSE, 
                               k = 400, gravity=1, ks=0.1, ksmax=10, delta = 1, center=NULL,
                               tolerate = 0.1, dim = 2, plotstep=10, plotlabels=TRUE){ 
  ####  g is a igraph network
  ####  iterations is the number of iterations to be performed
  ####  linlog is a variant which uses logarithmic attraction force F <- log (1+F)
  ####  pos is the table (NumberOfNodes x dimension) of the initial locations of points, if specified
  ####  nohubs is a variant in which nodes with high indegree have more central position than nodes with outdegree (for directed graphs) 
  ####  k is the repel constant : the greater the constant k the stronger the repulse force between points 
  ####  gravity is the gravity constant : indicates how strongly the nodes should be attracted to the center of gravity
  ####  ks is the speed constant : the greater the value of ks the more movement the nodes make under the acting forces
  ####  ksmax limits the speed from above
  ####  delta is the parameter to modify attraction force; means that weights are raised to the power = delta
  ####  center is the center of gravity  
  ####  tolerance is the tolerance to swinging constant
  ####  dim is the dimension
  ####  plotstep is the frequency of plotting intermediate iterations
  ####  plotlabels is TRUE if the labels should be included in the intermediate interations plot
  
  A <- igraph::get.adjacency(g, type="both",
                             attr=NULL, edges=FALSE, names=TRUE,
                             sparse=FALSE)
  
  #### center of gravity is by default set to the origin
  if(is.null(center)) center <- rep(0,dim)
  
  nnodes <- nrow(A)
  #### Binary will be a matrix of simple incidence (0-not connected, 1-connected)
  Binary <- A
  Binary[Binary!=0] <- 1
  #### Deg will be a vector of the degrees of vertices
  Deg <- rowSums(Binary)
  #### Forces1 will be a table containing all the sums of forces acting on points at a step
  Forces1 <- matrix(0, nrow = dim, ncol = nnodes)
  
  #### If there are no initial coordinates of points, 
  #### they are chosen at random from 1000^dim square
  if (is.null(pos))  {
    difference <- 2000/(nnodes*dim) 
    position <- matrix(sample(seq(-1000,1000,difference),nnodes*dim),nnodes,dim)
  }  else  {
    position <- pos
  }
  
  #### None of the nodes should be exactly at the center of gravity###
  temp <- which(position[,1] == center[1])
  for( index in 2:ncol(position)){
    temp <- intersect(temp,which(position[,index] == center[index]))
  }
  position[temp,] <- center + 0.01
  rm(index,temp)
  
  #### displacement will be a matrix of points' movement at the current iteration
  displacement <- matrix(rep(0,dim*nnodes),dim,nnodes)
  
  m <- nrow(position) 
  cur.step <- 0
  for (iteration in 1:iterations)
  {
    cur.step <- cur.step + 1
    cat(paste(Sys.time(), " step ", cur.step, "/", iterations, "\n", sep = ""))
    displacement <- displacement * 0
    #### Forces2 is the table of the forces from previous step
    #### Forces1 is the table of the forces from current step
    Forces2 <- Forces1 
    Forces1 <- matrix(, nrow = dim, ncol = 0)
    
    #### Calculate the Forces for each node
    ### Distance matrix between all nodes
    distances <- as.matrix(dist(position))
    distances[which(distances < 0.01)] <- 0.01 #We impose a minimum distance
    ### Each element of the list contains a matrix with the j = 1,2,..., dim dimension of the unitary vector 1    
    mylist <- vector("list",dim)
    for (j in 1:dim){
      mylist[[j]] <- (tcrossprod(position[,j],rep(1,m))-tcrossprod(rep(1,m),position[,j]))/distances  
    }
    ### Calculate the repulsion Force   
    Fr <- k*((tcrossprod(rep(1,m),Deg)+1)*(tcrossprod(Deg,rep(1,m))+1))/distances
    
    #The classical attraction force is just based on distance
    Fa <- distances   
    #The linlog mode calculates the attraction force as log(1+d(n1,n2))
    if(linlog){
      Fa <- log(1+Fa)
    }
    #Edge weights. The edges are weighted based on parameter delta. delta=0 implies no weight
    Fa <- (A^delta)*Fa
    
    #Dissuade Hubs. This mode is meant to grant authorities (nodes with high indegree)
    #a more central position than hubs (nodes with high outdegree)
    if(nohubs){
      Fa <- Fa/(tcrossprod(Deg,rep(1,m))+1)
    }
    
    ### Function to calculate the Attraction and Repulsion forces
    Farfunction <- function(x) rowSums(x*(Fr-Fa),na.rm=T)
    ### And we aggregate it over all dimensions
    Far <- do.call(rbind,lapply(mylist,Farfunction))
    ### Unitary Vector 2, the directions between each point and the center
    uv2 <- apply(matrix(rep(center,m),nrow=m,byrow=T)-position,1,function(x) x/sqrt(sum(x^2)))
    ### The gravity force
    #Fg <- uv2*matrix(rep(gravity*(rowSums(A)+1),dim),nrow=dim,byrow=T)
    Fg <- uv2*matrix(rep(gravity*(Deg+1),dim),nrow=dim,byrow=T)
    ### Forces 1 is the sum between all forces: Far (Fa + Fr) and Fg
    Forces1 <- Far+Fg
    Forces1 <- round(Forces1,2) #Use the first two decimals for the Forces.
    
    #### Swing is the vector of the swingings of all points
    swing <- abs(colSums((Forces1-Forces2)^2)^(1/2))
    Global_swing <- sum((Deg + 1)*swing)
    
    #### tra is the vector of the traction of all points
    tra <- abs(colSums((Forces1+Forces2)^2)^(1/2))/2
    Global_tra <- sum((Deg+1)*tra)
    
    #### Global speed calculation
    Global_speed <- tolerate * Global_tra/Global_swing
    #### speed is the vector of individual speeds of points
    speed <- ks * Global_speed /  (1 + Global_speed * (swing)^(1/2))
    
    #### Imposing constrains on speed
    speed_constrain <- ksmax/abs(colSums((Forces1^2))^(1/2))
    speed <- ifelse(speed>=speed_constrain,speed_constrain,speed)
    
    #### calculating displacement and final position of points after iteration
    displacement <- Forces1 * t(matrix(rep(speed,dim),nnodes,dim))
    position <- position + t(displacement)
    
    #### Iteration plot. This is simply to see the evolution of the positions over iterations
    #### Is much faster to visualize directly using R base plots instead of igraph plots
    
    if(!plotstep==0&dim==2){
      if(iteration%%plotstep==0)  {
        plot(position, main=paste0("iteration: ",iteration), xlab="", ylab="")
        if(plotlabels) text(position, labels=V(g)$name, cex= 0.7, pos=3)
      }
    }
  }
  return (position)  
}


