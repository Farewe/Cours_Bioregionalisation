scan_rmd_for_packages <- function(file_path) {
  # Lire le contenu du fichier Rmd
  lines <- readLines(file_path)
  
  # Rechercher les lignes qui contiennent library() ou require()
  library_lines <- grep("library\\(|require\\(", lines, value = TRUE)
  
  # Extraire les noms des packages en utilisant une expression régulière
  packages <- gsub(".*(library|require)\\((['\"]?)([^'\")]+)['\"]?\\).*", "\\3", library_lines)
  
  # Supprimer les doublons
  packages <- unique(packages)
  
  return(packages)
}

# Utiliser la fonction sur un fichier .Rmd spécifique
file_path <- c("correction2.Rmd",
               "fonctionc.Rmd")
packages <- unique(unlist(lapply(file_path, scan_rmd_for_packages)))

# Générer une ligne de commande unique pour installer les packages
install_command <- paste("install.packages(c(", paste0('"', packages, '"', collapse = ", "), "))")

# Afficher la ligne de commande d'installation
cat("Ligne de commande pour installer les packages :\n")
cat(install_command, "\n")
