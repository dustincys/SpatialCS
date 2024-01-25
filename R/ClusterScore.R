#' Calculate the cluster score for spatial transcriptomics data
#'
#' @importFrom readr read_tsv
#' @importFrom dplyr `%>%` tibble group_by count filter rename mutate
#' @param spatial_info_table
#' Path to a tab delimited text file with column: spot, x_array, y_array,
#' pathology_annotation, cluster, pathology_annotation_type or, a tibble object
#' that contains those columns
#' @param D distance
#' @param alpha weight paramter 1
#' @param beta weight paramter 2
#'
#' @return a score in range \[0,1\]
#' @export
#'
SpatialClusterScore <- function(spatial_info_table, D, alpha = 1, beta = 0) {
  if (is.character(spatial_info_table)) {
    sst <- read_tsv(spatial_info_table, show_col_types = F)
  } else {
    sst <- spatial_info_table
  }
  colnames(sst) <- c("spot", "x_array", "y_array", "pathology_annotation", "cluster")

  sst$q <- sst$pathology_annotation
  sst$q[sst$pathology_annotation == "unlabelled"] <- "Unlabelled"
  sst$q[is.na(sst$q)] <- "Unlabelled"

  contingency_matrix <- table(sst$q, sst$cluster)
  a <- a_fun(contingency_matrix[rownames(contingency_matrix)!="Unlabelled",])
  aplus1 <- a_fun(contingency_matrix)

  c <- c_fun(contingency_matrix[rownames(contingency_matrix)!="Unlabelled",], D, alpha, beta)

  if ((aplus1 + c) == 0 ) {
    return(1)
  } else {
    return(a / (aplus1 + c))
  }
}

a_fun <- function(cm) {
  if(is.matrix(cm)) {
    return(sum(apply(cm, MARGIN = c(1,2), FUN = function(x) {
      if(x>0){
        return(x*(x-1))
      }else{
        return(0)
      }
    })))
  } else {
    return(sum(sapply(cm, FUN = function(x) {
      if(x>0){
        return(x*(x-1))
      }else{
        return(0)
      }
    })))
  }
}

c_fun <- function(cm, D, alpha, beta) {
  c <- 0
  if (is.array(cm)) {
    for(i in 1:dim(cm)[1]) {
      for(j in 1:dim(cm)[1]) {
        i_name <- rownames(cm)[i]
        j_name <- rownames(cm)[j]
        c <- c + sum(cm[i,] * cm[j,] * (alpha * (D[i_name,j_name] + beta)))
      }
    }
  } else {
    for(i in 1:length(cm)) {
      for(j in 1:length(cm)) {
        i_name <- names(cm)[i]
        j_name <- names(cm)[j]
        c <- c + sum(cm[i] * cm[j] * (alpha * (D[i_name,j_name] + beta)))
      }
    }
  }
  return(c)
}

#' Calculate the distance between histology annotation
#'
#' @importFrom Seurat Load10X_Spatial SCTransform FindSpatiallyVariableFeatures SpatiallyVariableFeatures Idents AverageExpression Cells
#' @importFrom dplyr `%>%` tibble group_by count filter rename mutate
#' @importFrom readr read_tsv
#' @importFrom stats cor
#' @param cellranger_out
#' 10x spaceranger count folder
#' @param annotation_table
#' table contains pathology_annotation, cluster, spot columns
#' @return distance matrix
#' @export
get_spatial_corr <- function(cellranger_out, annotation_table){
  seurat_obj <- Seurat::Load10X_Spatial(cellranger_out)
  sst <- readr::read_tsv(annotation_table)
  colnames(sst) <- c("spot", "x_array", "y_array", "pathology_annotation", "cluster")
  seurat_obj <- subset(seurat_obj, cells = sst$spot)
  seurat_obj <- Seurat::SCTransform(seurat_obj, assay = "Spatial")
  seurat_obj <- Seurat::FindSpatiallyVariableFeatures(seurat_obj, assay = "SCT", selection.method = "moransi")
  svf <- SpatiallyVariableFeatures(seurat_obj, selection.method = "moransi")
  seurat_obj_svf <- subset(seurat_obj, features = svf)
  seurat_obj_svf$pathology_annotation <-
    sst$pathology_annotation[match(Cells(seurat_obj_svf), sst$spot)]

  Seurat::Idents(seurat_obj_svf) <- seurat_obj_svf$pathology_annotation
  avg.exp <- Seurat::AverageExpression(seurat_obj_svf, slot = "data")$Spatial
  avg.exp <- avg.exp[!is.na(colSums(avg.exp)),]
  cor.exp <- as.matrix(cor(avg.exp, method = "spearman"))

  # min_v <- min(cor.exp)
  # max_v <- max(cor.exp)
  # D <- (max_v - cor.exp) / (max_v - min_v)
  print(cor.exp)
  D <- 1 - cor.exp

  return(D)
}


#' Calculate the rogue score for spatial transcriptomics data
#'
#' @importFrom readr read_tsv read_csv
#' @importFrom Seurat Read10X_h5
#' @importFrom stringr str_replace
#' @importFrom ROGUE matr.filter SE_fun CalculateRogue
#'
#' @param exp_table
#' Path to a csv expression matrix file
#' @param spatial_info_table
#' Path to a tab delimited text file with column: spot, x_array, y_array,
#' pathology_annotation, cluster, pathology_annotation_type or, a tibble object
#' that contains those columns
#' @param min.cells minimum cell (spot) number to calculate ROGUE, default 10
#' @param min.genes minimum gene number to calculate ROGUE, default 10
#'
#' @return a score in range \[0,1\]
#' @export
#'
SpatialRogueScore <- function(exp_table, spatial_info_table, min.cells = 10, min.genes = 10) {
  if (is.character(spatial_info_table)) {
    sst <- read_tsv(spatial_info_table)
  } else {
    sst <- spatial_info_table
  }
  colnames(sst) <- c("spot", "x_array", "y_array", "pathology_annotation", "cluster")

  if (is.character(exp_table)) {
    exp <- as.matrix(Read10X_h5(exp_table, use.names = TRUE, unique.features = TRUE))
  } else {
    exp <- exp_table
  }

  pt <- table(sst %>%
               filter(! pathology_annotation %in%
                      c("Unlabelled", "unlabelled")) %>%
               pull(pathology_annotation))
  if(min(pt) < min.cells) {
    min.cells <- min(pt)
  }

  sst$rogue = 0
  for(ct in as.vector(unique(sst$cluster))){
    if(sum(sst$cluster == ct) < min.cells) {
      next
    }
    targetMatrix <- exp[, sst$spot[sst$cluster == ct]]
    if(! is.matrix(targetMatrix)) {
      next
    }
    if(dim(targetMatrix)[2] > 10000){
      targetMatrix <- targetMatrix[, sample(colnames(targetMatrix), 10000, replace = F)]
    }
    targetMatrix <- matr.filter(as.matrix(targetMatrix), min.cells = min.cells, min.genes = min.genes)
    if (dim(targetMatrix)[1] == 0) {
      next
    }
    targetSE <- SE_fun(targetMatrix)
    ct.rogue.value <- CalculateRogue(targetSE, platform = "UMI")
    sst$rogue[sst$cluster == ct] = ct.rogue.value
  }

  return(mean(sst$rogue))
}

#' Calculate the rogue score for spatial transcriptomics data
#'
#' @importFrom readr read_tsv read_csv
#' @importFrom Seurat Read10X_h5
#' @importFrom stringr str_replace
#' @importFrom ROGUE matr.filter SE_fun CalculateRogue
#'
#' @param exp_table
#' Path to a csv expression matrix file
#' @param spatial_info_table
#' Path to a tab delimited text file with column: spot, x_array, y_array,
#' pathology_annotation, cluster, pathology_annotation_type or, a tibble object
#' that contains those columns
#' @param min.cells minimum cell (spot) number to calculate ROGUE, default 10
#' @param min.genes minimum gene number to calculate ROGUE, default 10
#'
#' @return a score in range \[0,1\]
#' @export
#'
SpatialRogueScore_ha <- function(exp_table, spatial_info_table, min.cells = 10, min.genes = 10) {
  if (is.character(spatial_info_table)) {
    sst <- read_tsv(spatial_info_table)
  } else {
    sst <- spatial_info_table
  }
  colnames(sst) <- c("spot", "x_array", "y_array", "pathology_annotation", "cluster")

  if (is.character(exp_table)) {
    exp <- as.matrix(Read10X_h5(exp_table, use.names = TRUE, unique.features = TRUE))
  } else {
    exp <- exp_table
  }

  ha_s <- sst %>%
    filter(! pathology_annotation %in%
             c("Unlabelled", "unlabelled")) %>%
    pull(pathology_annotation) %>%
    unique()

  sst$ha_rogue = 0
  for(ha in ha_s){
    if(sum(sst$pathology_annotation == ha) < min.cells) {
      next
    }
    targetMatrix <- exp[, sst$spot[sst$pathology_annotation == ha]]
    if(! is.matrix(targetMatrix)) {
      next
    }
    if(dim(targetMatrix)[2] > 10000){
      targetMatrix <- targetMatrix[, sample(colnames(targetMatrix), 10000, replace = F)]
    }
    targetMatrix <- matr.filter(as.matrix(targetMatrix), min.cells = min.cells, min.genes = min.genes)
    if (dim(targetMatrix)[1] == 0) {
      next
    }
    targetSE <- SE_fun(targetMatrix)
    ha.rogue.value <- CalculateRogue(targetSE, platform = "UMI")
    sst$ha_rogue[sst$pathology_annotation == ha] = ha.rogue.value
  }

  # rt <- sst %>%
  #   filter(pathology_annotation %in% ha_s) %>%
  #   group_by(pathology_annotation) %>%
  #   summarize(mean_rogue = mean(ha_rogue))
  # rt <- rt %>%
  #   mutate(mean_rogue_all = mean(sst$ha_rogue))

  return(mean(sst$ha_rogue))
}


#' Get the plot ha cluster path to the Python script file in the package, and call it
#'
#' @param spatial_info_table cluster table path
#' @param color_table color table path
#' @param highlight_PAs_fill clusters to highlight in fill color
#' @param highlight_PAs_border clusters to highlight in border color
#' @param output_folder output folder path
#' @param output_file output file
#'
#' @export
plot_ha_cluster <- function(spatial_info_table,
                            color_table,
                            highlight_PAs_fill=NULL,
                            highlight_PAs_border=NULL,
                            output_folder = ".",
                            output_file="space_ha_c.pdf") {
  script_path <- system.file("python",
                             "plot_ha_cluster.py",
                             package="SpatialClusterScore")
  print(paste("script_path", script_path))


  if(is.null(highlight_PAs_fill)) {
    print(paste("python",
                script_path,
                paste0("--clusterTablePath=", paste0("\"", spatial_info_table, "\"")),
                paste0("--colorTablePath=", paste0("\"", color_table, "\"")),
                paste0("--outFileName=", paste0("\"", output_file, "\"")),
                paste0("--outFolderPath=", paste0("\"", output_folder, "\""))))
    system(paste("python",
                script_path,
                paste0("--clusterTablePath=", paste0("\"", spatial_info_table, "\"")),
                paste0("--colorTablePath=", paste0("\"", color_table, "\"")),
                paste0("--outFileName=", paste0("\"", output_file, "\"")),
                paste0("--outFolderPath=", paste0("\"", output_folder, "\""))))
  } else {


    print(paste("python",
                script_path,
                paste0("--clusterTablePath=", paste0("\"", spatial_info_table, "\"")),
                paste0("--colorTablePath=", paste0("\"", color_table, "\"")),
                paste0("--outFolderPath=", paste0("\"", output_folder, "\"")),
                paste0("--outFileName=", paste0("\"", output_file, "\"")),
                paste("--highlightFill",
                      paste0("\"", highlight_PAs_fill, "\"", collapse = " ")),
                paste("--highlightBorder",
                      paste0( "\"", highlight_PAs_border,  "\"", collapse = " "))
                ))
    system(paste("python",
                 script_path,
                 paste0("--clusterTablePath=", paste0("\"", spatial_info_table, "\"")),
                 paste0("--colorTablePath=", paste0("\"", color_table, "\"")),
                 paste0("--outFolderPath=", paste0("\"", output_folder, "\"")),
                 paste0("--outFileName=", paste0("\"", output_file, "\"")),
                 paste("--highlightFill",
                      paste0("\"", highlight_PAs_fill, "\"", collapse = " ")),
                 paste("--highlightBorder",
                      paste0( "\"", highlight_PAs_border,  "\"", collapse = " "))
                 ))
  }
}


#' Get the plot c involved path to the Python script file in the package, and call it
#'
#' @param spatial_info_table cluster table path
#' @param color_table color table path
#' @param exp_table
#' Path to a csv expression matrix file
#' @importFrom readr read_tsv
#' @importFrom Seurat Read10X_h5
#' @param highlight_PAs_fill clusters to highlight in fill color
#' @param highlight_PAs_border clusters to highlight in border color
#' @param output_folder output folder path
#' @param output_file output file
#' @param min_frac minimum fraction
#' @param min_n minimum number
#'
#' @export
plot_c_involved <- function(spatial_info_table,
                            color_table,
                            exp_table,
                            highlight_PAs_fill=NULL,
                            highlight_PAs_border=NULL,
                            min_frac = 0.1,
                            min_n = 10,
                            output_folder = ".",
                            output_file="space_ha_c_involved.pdf") {
  script_path <- system.file("python",
                             "plot_c_involved.py",
                             package="SpatialClusterScore")
  print(paste("script_path", script_path))


  if (is.character(spatial_info_table)) {
    sst <- read_tsv(spatial_info_table)
  } else {
    sst <- spatial_info_table
  }
  colnames(sst) <- c("spot", "x_array", "y_array", "pathology_annotation", "cluster")

  if (is.character(exp_table)) {
    exp <- as.matrix(Read10X_h5(exp_table, use.names = TRUE, unique.features = TRUE))
  } else {
    exp <- exp_table
  }

  c_involved <- sst %>%
    filter(pathology_annotation %in% highlight_PAs_fill) %>%
    group_by(cluster) %>%
    count() %>%
    ungroup() %>%
    mutate(frac = n / sum(n)) %>%
    filter(n > min_n | frac > min_frac) %>%
    pull(cluster) %>%
    unique()

  ha_average_rogue <- SpatialRogueScore(exp_table,
                                        sst %>% filter(pathology_annotation %in% highlight_PAs_fill),
                                        min.cells = 10, min.genes = 10)
  cluster_average_rogue <- SpatialRogueScore(exp_table,
                                             sst %>% filter(cluster %in% c_involved),
                                             min.cells = 10, min.genes = 10)

  if(is.null(highlight_PAs_fill)) {
    print(paste("python",
                script_path,
                paste0("--clusterTablePath=", paste0("\"", spatial_info_table, "\"")),
                paste0("--colorTablePath=", paste0("\"", color_table, "\"")),
                paste0("--outFileName=", paste0("\"", output_file, "\"")),
                paste0("--outFolderPath=", paste0("\"", output_folder, "\""))))
    system(paste("python",
                script_path,
                paste0("--clusterTablePath=", paste0("\"", spatial_info_table, "\"")),
                paste0("--colorTablePath=", paste0("\"", color_table, "\"")),
                paste0("--outFileName=", paste0("\"", output_file, "\"")),
                paste0("--outFolderPath=", paste0("\"", output_folder, "\""))))
  } else {
    print(paste("python",
                script_path,
                paste0("--haAvgRogue=", ha_average_rogue),
                paste0("--clusterAvgRogue=", cluster_average_rogue),
                paste0("--clusterTablePath=", paste0("\"", spatial_info_table, "\"")),
                paste0("--colorTablePath=", paste0("\"", color_table, "\"")),
                paste0("--outFolderPath=", paste0("\"", output_folder, "\"")),
                paste0("--outFileName=", paste0("\"", output_file, "\"")),
                paste0("--minFrac=", min_frac),
                paste0("--minN=", min_n),
                paste("--highlightFill",
                      paste0("\"", highlight_PAs_fill, "\"", collapse = " ")),
                paste("--highlightBorder",
                      paste0( "\"", highlight_PAs_border,  "\"", collapse = " "))
                ))
    system(paste("python",
                 script_path,
                 paste0("--haAvgRogue=", ha_average_rogue),
                 paste0("--clusterAvgRogue=", cluster_average_rogue),
                 paste0("--clusterTablePath=", paste0("\"", spatial_info_table, "\"")),
                 paste0("--colorTablePath=", paste0("\"", color_table, "\"")),
                 paste0("--outFolderPath=", paste0("\"", output_folder, "\"")),
                 paste0("--outFileName=", paste0("\"", output_file, "\"")),
                 paste0("--minFrac=", min_frac),
                 paste0("--minN=", min_n),
                 paste("--highlightFill",
                      paste0("\"", highlight_PAs_fill, "\"", collapse = " ")),
                 paste("--highlightBorder",
                      paste0( "\"", highlight_PAs_border,  "\"", collapse = " "))
                 ))
  }
}

