#' identify confident cells in each cell cluster
#'
#' @param object.list a list of Seurat object
#' @param features.integration features used to identify confident cells
#' @param quantile.cutoff quantile cutoff (default = 0.75)
#' @param min.cluster the minimum number of cells in the cluster for identifying confident cells
#' @param assay Assay to use
#' @importFrom stats quantile
#' @import Seurat
#'
#' @return a list of Seurat object with identified confident cells in each dataset
#' @export

identifyConfidentCells <- function(object.list, features.integration, quantile.cutoff = 0.75, min.cluster = 20, assay = NULL) {
  for (ii in 1:length(object.list)) {
    object <- object.list[[ii]]
    if (is.null(assay)) {assay <- DefaultAssay(object = object)}
    data.norm <- GetAssayData(object = object, assay = assay, slot = "data")
    data.norm <- data.norm[features.integration, ]
    graph.name <- paste0(assay, "_", "snn")
    SNN <- object[[graph.name]]
    cluster <- Idents(object)
    cluster.uni <- as.character(unique(cluster))

    flag <- c()
    for (i in 1:length(cluster.uni)) {
      Index <- which(cluster == cluster.uni[i])
      if (length(Index) > min.cluster) {
        Si <- SNN[Index, Index]
        colsi <- Matrix::colSums(Si)/nrow(Si)
        thresh <- stats::quantile(colsi,quantile.cutoff)
        index <- which(colsi > thresh)
        if (length(index) > min.cluster) {
          flag <- c(flag, Index[index])
        }
      }
    }
    confident_index <- flag
    confident.cells <- colnames(data.norm)[confident_index]
    Misc(object, slot = 'confident.cells') <- confident.cells
    object.list[[ii]] <- object
  }
  return(object.list)
}


#' Learn technical variation between any two datasets
#'
#' @param object.list List of Seurat objects
#' @param features.integration feature
#' @param similarity.cutoff a thresholding parameter determining whether cell clusters are shared across different datasets based on their similarity. If T is too small, the biological variation may be removed. If T is too large, the technical variation could not be removed.
#' @param assay assay to use
#' @import Seurat
#'
#' @return a list of two structured matrices: one is the data matrix and the other is the confounding matrix
#' @export
#'
#' @examples
learnTechnicalVariation <- function(object.list, features.integration, similarity.cutoff = 0.6, assay = NULL) {
  if (is.null(assay)) {assay <- DefaultAssay(object = object.list[[1]])}
  Num <- buildBipartiteGraph(object.list, features.integration, assay = assay)
  n <- length(object.list)
  idys <- list()
  c_record <- list()
  Vec <- c()
  flags <- c()
  for (i in 1:length(object.list)) {
    object <- object.list[[i]]
    confident.cells <-  Misc(object, slot = 'confident.cells')
    subcluster <- Idents(object)[confident.cells]
    c <- sort(unique(subcluster)) # order
    idy <- list()
    for (j in 1:length(c))
    {
      idy[[j]] <- which(subcluster == c[[j]])
    }
    idys[[i]] <- idy
    c_record[[i]] <- seq(from = 1, to = length(c), by = 1)
    Vec[i] <- length(c)
    flags <- c(flags, matrix(i, nrow = 1, ncol = length(c)))
  }

  S <- cumsum(Vec)
  mm <- which(Num > similarity.cutoff, arr.ind = TRUE)
  idy1 <- mm[,1]
  idy2 <- mm[,2]
  # If one column has more than one corresponding rows, take both from different conditions, each condition take the largest value
  hh <- table(idy1)
  m1 <-  as.numeric(hh)
  n1 <- sort(unique(idy1))
  n1_group <- list()
  n1_group[[1]] <- n1[which(n1 <= S[1])]
  m1_group <- list()
  m1_group[[1]] <- m1[which(n1 <= S[1])]

  idy1_group <- list()
  idy1_group[[1]] <- idy1[which(idy1 <= S[1])]
  idy2_group <- list()
  idy2_group[[1]] <- idy2[which(idy1 <= S[1])]

  for (i in 2:length(object.list)) {
    id1 <- intersect(which(n1 > S[i-1]), which(n1 <= S[i]))
    n1_group[[i]] <- n1[id1]
    m1_group[[i]] <- m1[id1]

    id2 <- intersect(which(idy1 > S[i-1]), which(idy1 <= S[i]))
    idy1_group[[i]] <- idy1[id2]
    idy2_group[[i]] <- idy2[id2]
  }

  # record index to form Ys and positions
  Idy1_record <- list()
  Idy2_record <- list()
  for (i in 1:n) {
    idy1_new <- c()
    idy2_new <-c()
    m1 <- m1_group[[i]]
    n1 <- n1_group[[i]]
    idy1 <- idy1_group[[i]]
    idy2 <- idy2_group[[i]]
    for (p in 1:length(m1)) {
      if (m1[p] > 1) {
        p1 <- which(idy1 == n1[p])
        # see whether they are in the same group
        candidate <- idy2[p1]
        flag_i <- flags[candidate] # group information
        # divided groups, each group select largest one
        x <- sort(unique(flag_i)) # group information
        for (q in 1:length(x)) {
          pq <- which(flag_i == x[q])# group index
          # return index
          idy2_s <- idy2[p1]
          idy2_s <- idy2_s[pq]
          values <- Num[n1[p],idy2_s] # take out values
          q1 <- which(values == max(values)) # return the largest index
          idy1_new <- c(idy1_new,idy1[p1[q1[1]]])
          idy2_new <- c(idy2_new,idy2_s[q1[1]])
        }
      } else {
        p2 <- which(idy1 == n1[p])
        idy1_new <- c(idy1_new,idy1[p2])
        idy2_new <- c(idy2_new,idy2[p2])
      }
    }
    # divided into different groups
    idy1_record <- list()
    idy2_record <- list()
    for (j in (i+1):n) {
      x1 <- which(flags[idy2_new] == j)
      idy1_record[[j]] <- idy1_new[x1]
      idy2_record[[j]] <- idy2_new[x1]
    }
    Idy1_record[[i]] <- idy1_record
    Idy2_record[[i]] <- idy2_record
  }
  # build X and Ys
  # X according to each cluster
  X_new_all <- matrix(0, nrow = length(features.integration), ncol =  1)
  vecN_all <- c()
  indexX_all <- c()
  for (j in 1:n) {
    X_new <- matrix(0, nrow = length(features.integration), ncol =  1)
    vecN <- c()
    indexX <- c()
    c <- c_record[[j]]
    IDx <- idys[[j]]
    confident.cells <- Misc(object.list[[j]], slot = 'confident.cells')
    data.use <- GetAssayData(object = object.list[[j]], assay = assay, slot = "data")
    X <- data.use[features.integration, confident.cells, drop = FALSE]

    for (k in 1:length(IDx)) {
      id <- which(c == k)
      vecN[k] <- length(IDx[[id]])
      X_new <- cbind(X_new, X[,IDx[[id]]])
      indexX <- c(indexX,IDx[[id]])
    }
    X_new <- X_new[,-1]
    X_new_all <- cbind(X_new_all, X_new)
    vecN_all <- c(vecN_all, vecN)
    indexX_all <- c(indexX_all, indexX)
  }

  X_b <- X_new_all[,-1]
  Sum_all <- cumsum(vecN_all)
  N <- Sum_all[length(vecN_all)]
  # build Y, Yi according to corresponding relationship
  Y_b <- matrix(0, nrow = ncol(X_b), ncol = 1)
  for (p in 1:(n-1)) {
    Y_new_p <- list()
    for (q in (p+1):n) {
      Idy1_p <- Idy1_record[[p]]
      Idy2_p <- Idy2_record[[p]]
      Idy1 <- Idy1_p[[q]]
      Idy2 <- Idy2_p[[q]]
      Ys <- list()
      if (length(Idy1) > 0) {
        for (i in 1:length(Idy1)) {
          val <- min(vecN_all[Idy1[i]], vecN_all[Idy2[i]])
          y1 <- diag(val)
          y2 <- -1*diag(val)
          if (Idy1[i] == 1) {
            Y1 <- rbind(y1, matrix(0, nrow = N-nrow(y1), ncol = val))
          } else {
            ya <- matrix(0, nrow = Sum_all[Idy1[i]-1], ncol = val)
            yb <- matrix(0, nrow = N-Sum_all[Idy1[i]-1]-nrow(y1), ncol = val)
            Y1 = rbind(ya, rbind(y1, yb))
          }
          yc <- matrix(0, nrow = Sum_all[Idy2[i]-1], ncol = val)
          yd <- matrix(0,nrow = N-Sum_all[Idy2[i]-1]-nrow(y2), ncol = val)
          Y2 = rbind(yc, rbind(y2, yd))
          Ys[[i]] <- Y1+Y2
        }
        Y_new <- matrix(0, nrow = nrow(Y1), ncol = 1)
        for (j in 1:length(Idy1)) {
          Y_new <- cbind(Y_new,Ys[[j]])
        }
        Y_new <- Y_new[,-1]
        Y_b <- cbind(Y_b, Y_new)
      }
    }
  }
  Y_b <- Y_b[,-1]
  return(list(X = X_b, Y = Y_b))
}



#' build Bipartite Graph
#'
#' @param object.list List of Seurat objects
#' @param features.integration feature to use
#' @param assay assay to use
#' @import Seurat
#'
#' @return a Bipartite Graph
buildBipartiteGraph <- function(object.list, features.integration, assay = NULL) {
  idys <- list()
  gene_groups <- list()
  t <- 0
  vec <- c()
  subXs <- c()
  for (i in 1:length(object.list)){
    object <- object.list[[i]]
    confident.cells <- Misc(object, slot = 'confident.cells')
    subcluster <- Idents(object)[confident.cells]
    data.use <- GetAssayData(object = object, assay = assay, slot = "data")
    subX <- data.use[features.integration, confident.cells, drop = FALSE]
    subXs <- cbind(subXs, subX)

    c <- as.character(sort(unique(subcluster))) # order
    subcluster <- as.character(subcluster)
    idy <- list()
    for (j in 1:length(c)) {
      idy[[j]] <- which(subcluster == c[j])+t
    }
    idys <- c(idys,idy)
    t <- t + ncol(subX)
    # marker
    markers <-  Misc(object, slot = 'markers')
    if (!is.null(markers)){
      identity <- as.character(markers$cluster)
      gene <- markers$gene
      gene_group <- list()
      for (j in 1:length(unique(subcluster))){
        if (sum(identity == c[j]) > 0) {
          gene_group[[j]] <- gene[identity == c[j]]# note: use int identity
        } else {
          s_subX <- Matrix::rowSums(subX[ ,subcluster == c[j]])
          z <- scale(s_subX)
          gene_group[[j]] <- rownames(subX)[z > 5]
        }
      }
      gene_groups <- c(gene_groups,gene_group)
      vec <- c(vec,length(unique(subcluster)))
    } else {
      gene_group <- list()
      for (j in 1:length(unique(subcluster))) {
        s_subX <- Matrix::rowSums(subX[,subcluster == c[j]])
        z <- scale(s_subX)
        gene_group[[j]] <- rownames(subX)[z > 5]
      }
      gene_groups <- c(gene_groups,gene_group)
      vec <- c(vec,length(unique(subcluster)))
    }
  }


  # varied jaccard
  S <- sum(vec)
  Num1 <- matrix(0,nrow = S,ncol = S)
  for (i in 1:(S-1)){
    for (j in (i+1):S){
      Num1[i,j] <- length(intersect(gene_groups[[i]],gene_groups[[j]]))/(min(length(gene_groups[[i]]),length(gene_groups[[j]])))
    }
  }
  Num1 <- Num1+t(Num1)
  Num1[is.na(Num1)] <- 0
  # cosine
  Num2 <- matrix(0,nrow = S,ncol = S)
  for (i in 1:(S-1)){
    for (j in (i+1):S){
      sxi <- subXs[,idys[[i]]]
      sxj <- subXs[,idys[[j]]]
      msxi <- as.matrix(apply(sxi,1,mean))# row mean
      msxj <- as.matrix(apply(sxj,1,mean))
      Num2[i,j] <- sum(msxi*msxj)/(norm(msxi,type = "F")%*%norm(msxj, type = "F"))
    }
  }
  Num2 <- Num2+t(Num2)
  # average
  Num = (Num1+Num2)/2
  return(Num)
}



#' Infer a shared embedding of cells across all datasets after removing technical variation
#'
#' @param object unintegrated Seurat object
#' @param mat structured matrices obtained from learnTechnicalVariation
#' @param assay assay to store the integrated data
#' @param new.assay.name Name for the new assay containing the integrated data
#' @param nDims.scMC number of dimensions to compute in the scMC integrated space
#' @param lambda the tuning parameter, non-negative.
#' @importFrom RSpectra eigs_sym
#' @import Seurat

#' @return Seurat object with an integrated assay
#' @export
#'
#'
integrateData <- function(object, mat, assay = NULL, new.assay.name = NULL, nDims.scMC = 40, lambda = 1){
  nDims <- nDims.scMC
  if (is.null(assay)) {assay <- DefaultAssay(object = object)}
  X <- t(as.matrix(mat$X))
  Y <- mat$Y
  nsam <- dim(X)[1]
  p <- dim(X)[2]
  if (is.null(dim(Y))){
    Y <- matrix(Y, ncol=1)
  }

  centerX = TRUE; centerY = TRUE; scaleX = FALSE; scaleY = FALSE;
  X <- scale(X, center = centerX, scale = scaleX)
  Y <- scale(Y, center = centerY, scale = scaleY)

  # calculate kernel matrix for Y
  K <- tcrossprod(Y)

  # learn the correction vectors by performing eigen-decomposition
  result_acpca <- RSpectra::eigs_sym(calAv, k=nDims, which = "LA", n=p, args=list(X=X, K=K, lambda=lambda))
  v <- matrix(result_acpca$vectors, ncol=nDims)
  data.joint <- GetAssayData(object = object, assay = assay, slot = "data")
  data.joint <- data.joint[VariableFeatures(object), ]
  emb <- crossprod(as.matrix(data.joint), v)
  emb <- as.matrix(emb)
  colnames(emb) <- paste0("scMC", seq_len(ncol(emb)))

  if (!is.null(new.assay.name)) {
    integrated.assay <- CreateAssayObject(data = data.joint)
    object[[new.assay.name]] <- integrated.assay
    assay <- new.assay.name
    DefaultAssay(object) <- assay
  }
  key <- paste0('scMC', '_')
  embeddings <- emb
  dr <- CreateDimReducObject(embeddings = embeddings, assay = assay, key = key)
  object[['scMC']] <- dr

  Misc(object, slot = 'correction.vectors') <- v
  Misc(object, slot = 'parameter')$lambda <- lambda

  return(object)
}

#' compute the projected data
#' @param v correction vectors
#' @param args arguments
#'
calAv <- function(v, args){
  X <- args$X
  K <- args$K
  lambda <- args$lambda
  return(crossprod(X, (diag(dim(K)[1])-lambda*K)%*%(X%*%matrix(v, ncol=1))) )
}


#######################################################
#######################################################
#' Identify cell clusters
#'
#' @param object a list of Seurat objects or a single Seurat object
#' @param mode "separate" or "integrated"; "separate" means performing clustering for each dataset separately, and "integrated" means performing clustering on the integrated data
#' @param resolution the resolution in Leiden algorithm; if it is NULL, the optimal resoultion will be inferred based on eigen spectrum
#' @param method Method for running leiden (defaults to matrix which is fast for small datasets). Enable method = "igraph" to avoid casting large data to a dense matrix.
#' @param algorithm Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm). Leiden requires the leidenalg python.
#' @param resRange the range of resolution values in Leiden algorithm; if it is NULL, the default range of resoultion will be from 0.1 to 0.5
#' @param nDims.consensus the number of singular values to estimate from the consensus matrix.
#' @param clustering.method method for performing clustering on the consensus matrix from a range of resolutions
#' @param graph.name Name of graph to use for the clustering algorithm
#' @param ... other parameter passing to FindClusters function in Seurat package
#' @importFrom Seurat Idents
#'
#' @return a list of Seurat objects or a single Seurat object
#' @export

identifyClusters <- function(object, mode = c("separate","integrated"), resolution = NULL,  method = c("matrix", "igraph"),algorithm = 4, resRange = NULL, nDims.consensus = 30, clustering.method = c("hierarchical", "community"), graph.name = NULL,...){
  mode <- match.arg(mode)
  method <- match.arg(method)
  clustering.method <- match.arg(clustering.method)

  if (mode == "separate") {
    for (i in 1:length(object)){
      cat("Identifying cell clusters in Dataset", i, '\n')
      object.i = object[[i]]
      if (is.null(graph.name)) {graph.name <- paste0(DefaultAssay(object = object.i), "_snn")}
      object.i <- identifyClusters_internal(object = object.i, resRange = resRange, resolution = resolution, method = method, nDims.consensus = nDims.consensus, clustering.method = clustering.method, graph.name = graph.name, ...)
      object[[i]] <- object.i
      cat("The initial guessed number of clusters is", nlevels(Idents(object.i)), '\n')
    }
  } else if (mode == "integrated") {
    cat("Identifying cell clusters in the integrated space", '\n')
    if (is.null(graph.name)) {graph.name <- paste0(DefaultAssay(object = object), "_snn")}
    object <- identifyClusters_internal(object = object, resRange = resRange, resolution = resolution, method = method, nDims.consensus = nDims.consensus, clustering.method = clustering.method, graph.name = graph.name,...)
    cat("The number of clusters in the integrated data is", nlevels(Idents(object)), '\n')
  }
  return(object)
}


#' identify cell clusters of a single Seurat object
#'
#' @param object Seurat object
#' @param resolution the resolution in Leiden algorithm; if it is NULL, the optimal resoultion will be inferred based on eigen spectrum
#' @param method Method for running leiden (defaults to matrix which is fast for small datasets). Enable method = "igraph" to avoid casting large data to a dense matrix.
#' @param algorithm Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm). Leiden requires the leidenalg python.
#' @param resRange the range of resolutions in Leiden algorithm; if it is NULL, the optimal range of resoultion will be from 0.1 to 0.5
#' @param nDims.consensus the number of singular values to estimate from the consensus matrix.
#' @param clustering.method method for performing clustering on the consensus matrix from a range of resolutions
#' @param graph.name Name of graph to use for the clustering algorithm
#' @param ... other parameter passing to FindClusters function in Seurat package
#' @importFrom  Seurat FindClusters Idents
#' @importFrom irlba irlba
#' @importFrom stats as.dist hclust cutree
#' @importFrom future plan nbrOfWorkers
#' @importFrom pbapply pbsapply
#' @importFrom future.apply future_sapply
#' @importFrom Matrix Matrix
#' @return Seurat object
#' @export
identifyClusters_internal <- function(object, resolution = NULL, method = "matrix", algorithm = 4, resRange = NULL, nDims.consensus = 30, clustering.method = c("hierarchical", "community"), graph.name = NULL,...) {
  clustering.method <- match.arg(clustering.method)
  if (is.null(graph.name)) {graph.name <- paste0(DefaultAssay(object = object), "_snn")}
  if (!is.null(resolution)) {
    object <- FindClusters(object, resolution = resolution, method = method, algorithm = algorithm, graph.name = graph.name, verbose = FALSE, ...)
  } else {
    N <- ncol(object)
    if (is.null(resRange)) {
      #resRange <- seq(0.1,0.5,by = 0.1)
      # resRange <- c(seq(0.1,1,by = 0.2), seq(1,2,by = 0.3))
      resRange <- c(c(0.02, 0.06), seq(0.1,1.5,by = 0.2))
    }

    object <- FindClusters(object, resolution = resRange, method = method, algorithm = algorithm, graph.name = graph.name, verbose = FALSE, ...)
    clustering_results <- object@meta.data[, paste0(graph.name, "_res.", resRange)]
    CM <- Matrix::Matrix(0, nrow = N, ncol = N, sparse = TRUE)
    ncluster <- c()
    for (i in 1:length(resRange)) {
      idents <- clustering_results[, i]
      clusIndex <- as.numeric(as.character(idents))
      CM <- CM + Matrix::Matrix(as.numeric(outer(clusIndex, clusIndex, FUN = "==")), nrow = N, ncol = N)
      ncluster[i] <- length(unique(idents))
    }
    CM <- CM/length(resRange)
    nCell = nrow(CM)
    num = length(unique(CM@x))
    if (num > 1) {
      # dertermine K
      nPC <- min(nDims.consensus, nCell)
      out <- irlba::irlba(CM, nv = nPC)
      s <- out$d
      # compute ratio
      cs = cumsum(s)/sum(s)
      tol <- 0.01
      K <- min(which((s/sum(s) < tol) == 1))-1

      if (clustering.method == "hierarchical") {
        cat("Peforming hierarchical clustering on the computed consensus matrix \n")
        d <- stats::as.dist(1 - CM)
        # hc <- stats::hclust(d, "ave")
        hc <- fastcluster::hclust(d, "ave")
        idents<- stats::cutree(hc,k = K)
        idents <- as.factor(idents)
        Idents(object) <- idents
      } else if (clustering.method == "community")  {
        resolution <- resRange[min(which(abs(ncluster - K) == min(abs(ncluster - K))))]
        cat("Peforming clustering using community-dection method with a resolution", resolution, "\n")
        object <- FindClusters(object, resolution = resolution, method = method, algorithm = algorithm, verbose = FALSE,...)
      }
    }
    else {
      idents <- as.factor(as.vector(matrix(1, nrow = nCell, ncol = 1)))
      Idents(object) <- idents
    }
  }
  levels <- levels(x = object)
  levels <- tryCatch(
    expr = as.numeric(x = levels),
    warning = function(...) {
      return(levels)
    },
    error = function(...) {
      return(levels)
    }
  )

  Idents(object = object) <- factor(x = Idents(object = object), levels = sort(x = levels))
  return(object)
}


#' identify shared nearest neighbors
#'
#' @param object a lsit of Seurat objects or a single Seurat object
#' @param mode "separate" or "integrate"
#' @param nDims.pca the number of dimensions to use for running PCA
#' @param force.pca Set force.pca = FALSE to skip the PCA calculation. Default = TRUE will calculate PCA.
#' @param nDims.knn the number of dimensions to use for building SNN
#' @param k.param Defines k for the k-nearest neighbor algorithm
#' @param prune.SNN Sets the cutoff for acceptable Jaccard index when computing the neighborhood overlap for the SNN construction. Any edges with values less than or equal to this will be set to 0 and removed from the SNN graph. Essentially sets the strigency of pruning (0 — no pruning, 1 — prune everything).
#' @param ... other parameters in FindNeighbors
#' @importFrom Seurat RunPCA FindNeighbors
#'
#' @return seurat object
#' @export

identifyNeighbors <- function(object, mode = c("separate","integrated"), nDims.pca = 40, force.pca = TRUE, nDims.knn = 40, k.param = 20, prune.SNN = 1/15, ...){
  mode <- match.arg(mode)
  if (mode == "separate") {
    for (i in 1:length(object)){
      object.i = object[[i]]
      if (!force.pca & ("pca" %in% names(object.i@reductions))) {
        message("Skiping the PCA calculation due to the existing PCA space")
      } else {
        cat("Performing PCA in Dataset ", i, '\n')
        if ("pca" %in% names(object.i@reductions)) {
          message("An existing PCA space is detected! Set force.pca = FALSE to skip the PCA calculation")
        }
        object.i <- RunPCA(object.i, npcs = nDims.pca, verbose = FALSE)
      }
      cat("Building SNN on the PCA space in Dataset ", i, '\n')
      object.i <- FindNeighbors(object.i, reduction = "pca", dims = 1:min(nDims.knn, nDims.pca), k.param = k.param, prune.SNN = prune.SNN, verbose = FALSE, ...)
      object[[i]] <- object.i
    }
  } else if (mode == "integrated") {
    cat("Building SNN on the integrated space", '\n')
    object <- FindNeighbors(object, reduction = "scMC", dims = 1:min(nDims.knn, nDims.pca),k.param = k.param, prune.SNN = prune.SNN, verbose = FALSE, ...)
  }
  return(object)
}

#' Identify marker genes associated with cell clusters
#'
#' @param object a lsit of Seurat objects or a single Seurat object
#' @param mode "separate" or "integrated"; "separate" means identifying marker genes of cell clusters from each dataset separately, and "integrated" means identifying marker genes of cell clusters from the integrated data
#' @param features features used to perform statistical test
#' @param test.use which test to use
#' @param only.pos Only return positive markers
#' @param min.pct Threshold of the percent of cells enriched in one cluster
#' @param logfc.threshold Threshold of Log Fold Change
#' @param ... other parameters in FindAllMarkers
#'
#' @import Seurat
#' @export
#'
#' @examples
identifyMarkers <- function(object, mode = c("separate","integrated"), features = NULL, test.use = "wilcox",
                            only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, ...) {
  mode <- match.arg(mode)
  if (mode == "separate") {
    for (i in 1:length(object)){
      object.i <- object[[i]]
      confident.cells <-  Misc(object.i, slot = 'confident.cells')
      object.subset <- subset(object.i, cells = confident.cells)
      cat("Identifying marker genes of cell clusters in Dataset", i, '\n')
      markers <- FindAllMarkers(object.subset, features = features, only.pos = only.pos, min.pct = min.pct, logfc.threshold = logfc.threshold, ...)
      Misc(object.i, slot = 'markers') <- markers
      object[[i]] <- object.i
    }
  } else if (mode == "integrated") {
    markers <- FindAllMarkers(object, features = features, only.pos = only.pos, min.pct = min.pct, logfc.threshold = logfc.threshold, ...)
    Misc(object, slot = 'markers') <- markers
  }
  return(object)
}

#' Identify integration features
#'
#' @param object.list a lsit of Seurat objects
#' @param integrationFeatures.method "joint" or "individual";
#' "joint": Identify integration features from the concatenated data matrix;
#' "individual": ranks features by the number of datasets they are deemed variable in, breaking ties by the median variable feature rank across datasets. It returns
#' the top scoring features by this ranking.
#' @param selection.method The method to choose top variable features:
#' vst: First, fits a line to the relationship of log(variance) and log(mean) using local polynomial regression (loess). Then standardizes the feature values using the observed mean and expected variance (given by the fitted line). Feature variance is then calculated on the standardized values after clipping to a maximum (see clip.max parameter).
#' mean.var.plot (mvp): First, uses a function to calculate average expression (mean.function) and dispersion (dispersion.function) for each feature. Next, divides features into num.bin (deafult 20) bins based on their average expression, and calculates z-scores for dispersion within each bin. The purpose of this is to identify variable features while controlling for the strong relationship between variability and average expression.
#' @param nfeatures Number of features to select as top variable features; only used when selection.method is set to 'vst'
#' @param mean.cutoff A two-length numeric vector with low- and high-cutoffs for feature means
#' @param dispersion.cutoff A two-length numeric vector with low- and high-cutoffs for feature dispersions
#' @param ... other parameters in FindVariableFeatures or SelectIntegrationFeatures
#' @import Seurat
#'
#' @export
#' @return a char vector containing the features for integration analysis
#'
identifyIntegrationFeatures <- function(object.list, integrationFeatures.method = c("joint", "individual"), selection.method = c("vst", "mean.var.plot"), nfeatures = 2000, mean.cutoff = c(0.01, 5), dispersion.cutoff = c(0.25, Inf), ...) {
  integrationFeatures.method <- match.arg(integrationFeatures.method)
  selection.method <- match.arg(selection.method)
  if (integrationFeatures.method == "joint") {
    options(warn = -1)
    object <- merge(x = object.list[[1]],y = object.list[2:length(x = object.list)])
    object <- FindVariableFeatures(object, selection.method = selection.method, nfeatures = nfeatures, mean.cutoff = mean.cutoff, dispersion.cutoff = dispersion.cutoff, verbose = FALSE, ...)
    features.integration <- VariableFeatures(object)
  } else (integrationFeatures.method == "individual")
  features.integration <- SelectIntegrationFeatures(object.list, nfeatures = nfeatures, ...)
}

#' Seurat wrapper for scMC
#'
#' Run scMC algorithm with Seurat pipelines
#'
#' @param object.list a list of Seurat objects, one per dataset
#' Parameters in identifyClusters
#' @param resolution the resolution in Leiden algorithm; if it is NULL, the optimal resoultion will be inferred based on eigen spectrum
#' @param method Method for running leiden (defaults to matrix which is fast for small datasets). Enable method = "igraph" to avoid casting large data to a dense matrix.
#' @param algorithm Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm). Leiden requires the leidenalg python.
#' @param resRange the range of resolution values in Leiden algorithm; if it is NULL, the default range of resoultion will be from 0.1 to 0.5
#' @param nDims.consensus the number of singular values to estimate from the consensus matrix.
#' @param clustering.method method for performing clustering on the consensus matrix from a range of resolutions
#' @param graph.name Name of graph to use for the clustering algorithm
#'
#' Parameters in identifyConfidentCells
#' @param quantile.cutoff quantile cutoff (default = 0.75)
#' Parameters in learnTechnicalVariation
#' @param similarity.cutoff a thresholding parameter determining whether cell clusters are shared across different datasets based on their similarity. If T is too small, the biological variation may be removed. If T is too large, the technical variation could not be removed.
#' Parameters in integrateData
#' @param new.assay.name Name for the new assay containing the integrated data
#' @param nDims.scMC number of dimensions to compute in the scMC integrated space
#' @param lambda the tuning parameter, non-negative.
#'
#' Parameters in identifyIntegrationFeatures
#' @param integrationFeatures.method "joint" or "individual";
#' "joint": Identify integration features from the concatenated data matrix;
#' "individual": ranks features by the number of datasets they are deemed variable in, breaking ties by the median variable feature rank across datasets. It returns
#' the top scoring features by this ranking.
#' @param selection.method The method to choose top variable features:
#' vst: First, fits a line to the relationship of log(variance) and log(mean) using local polynomial regression (loess). Then standardizes the feature values using the observed mean and expected variance (given by the fitted line). Feature variance is then calculated on the standardized values after clipping to a maximum (see clip.max parameter).
#' mean.var.plot (mvp): First, uses a function to calculate average expression (mean.function) and dispersion (dispersion.function) for each feature. Next, divides features into num.bin (deafult 20) bins based on their average expression, and calculates z-scores for dispersion within each bin. The purpose of this is to identify variable features while controlling for the strong relationship between variability and average expression.
#' @param nfeatures Number of features to select as top variable features; only used when selection.method is set to 'vst'
#' @param mean.cutoff A two-length numeric vector with low- and high-cutoffs for feature means
#' @param dispersion.cutoff A two-length numeric vector with low- and high-cutoffs for feature dispersions
#'
#' Parameters in identifyNeighbors
#' @param nDims.pca the number of dimensions to use for running PCA
#' @param force.pca Set force.pca = FALSE to skip the PCA calculation. Default = TRUE will calculate PCA.
#' @param nDims.knn the number of dimensions to use for building SNN
#' @param k.param Defines k for the k-nearest neighbor algorithm
#' @param prune.SNN Sets the cutoff for acceptable Jaccard index when computing the neighborhood overlap for the SNN construction. Any edges with values less than or equal to this will be set to 0 and removed from the SNN graph. Essentially sets the strigency of pruning (0 — no pruning, 1 — prune everything).
#'
#' Parameters in identifyMarkers
#' @param features features used to perform statistical test
#' @param test.use which test to use
#' @param only.pos Only return positive markers
#' @param min.pct Threshold of the percent of cells enriched in one cluster
#' @param logfc.threshold Threshold of Log Fold Change
#'
#' @param add.cell.ids A character vector of length(object.list) when merging multiple objects. Appends the corresponding values to the start of each objects' cell names.
#' @param assay Assay to use
#' @param ... other parameter passing to Seurat functions
#' @return A Seurat object with the integrated space from scMC
#' @export
RunscMC <- function(object.list,
                    resolution = NULL,  method = c("matrix", "igraph"), algorithm = 4, resRange = NULL, nDims.consensus = 30, clustering.method = c("hierarchical", "community"), graph.name = NULL,# parameters in identifyClusters
                    quantile.cutoff = 0.75, # parameters in identifyConfidentCells
                    similarity.cutoff = 0.6, # parameters in learnTechnicalVariation
                    new.assay.name = NULL, nDims.scMC = 40, lambda = 1, # parameters in integrateData
                    integrationFeatures.method = c("joint", "individual"), selection.method = c("vst", "mean.var.plot"), nfeatures = 2000, mean.cutoff = c(0.01, 5), dispersion.cutoff = c(0.25, Inf), # parameters in identifyIntegrationFeatures
                    nDims.pca = 40, force.pca = TRUE, nDims.knn = 40, k.param = 20, prune.SNN = 1/15, # parameters in identifyNeighbors
                    features = NULL, test.use = "wilcox", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, # parameters in identifyMarkers
                    add.cell.ids = NULL, assay = "RNA",
                    ...) {
  method <- match.arg(method)
  integrationFeatures.method <- match.arg(integrationFeatures.method)
  selection.method <- match.arg(selection.method)
  clustering.method <- match.arg(clustering.method)
  if (is.null(graph.name)) {graph.name <- paste0(assay, "_snn")}
  # step2. identify clusters with different resolution for each condition
  # compute SNN
  object.list <- tryCatch({
    identifyNeighbors(object.list, nDims.pca = nDims.pca, nDims.knn = nDims.knn, k.param = k.param, prune.SNN = prune.SNN, ...)
  }, error = function(e) {
    identifyNeighbors(object.list, nDims.pca = nDims.pca, nDims.knn = nDims.knn, k.param = k.param, prune.SNN = prune.SNN)
  })
  # identify clusters
  object.list <- tryCatch({
    identifyClusters(object.list, resolution = resolution,  method = method, algorithm = algorithm, resRange = resRange, nDims.consensus = nDims.consensus, clustering.method = clustering.method, graph.name = graph.name, ...)
  }, error = function(e) {
    identifyClusters(object.list, resolution = resolution,  method = method, algorithm = algorithm, resRange = resRange, nDims.consensus = nDims.consensus, clustering.method = clustering.method, graph.name = graph.name)
  })

  ## step3. detect cluster-specific cells with high confident
  features.integration <- tryCatch({
    identifyIntegrationFeatures(object.list, selection.method = selection.method, nfeatures = nfeatures, mean.cutoff = mean.cutoff, dispersion.cutoff = dispersion.cutoff, ...)
  }, error = function(e) {
    identifyIntegrationFeatures(object.list, selection.method = selection.method, nfeatures = nfeatures, mean.cutoff = mean.cutoff, dispersion.cutoff = dispersion.cutoff)
  })
  object.list <- identifyConfidentCells(object.list, features.integration, quantile.cutoff = quantile.cutoff, assay = assay)

  ## step4. Identify marker genes associated with the putative cell clusters in each dataset
  object.list <- tryCatch({
    identifyMarkers(object.list, features = features, test.use = test.use,
                    only.pos = only.pos, min.pct = min.pct, logfc.threshold = logfc.threshold, ...)
  }, error = function(e) {
    identifyMarkers(object.list, features = features, test.use = test.use,
                    only.pos = only.pos, min.pct = min.pct, logfc.threshold = logfc.threshold)
  })

  ## step 5. Learn technical variation between any two datasets
  structured_mat <- learnTechnicalVariation(object.list, features.integration, similarity.cutoff = similarity.cutoff, assay = assay)

  ## step 7. Learn a shared embedding of cells across all datasets after removing technical variation
  object <- merge(x = object.list[[1]],y = object.list[2:length(x = object.list)], add.cell.ids = add.cell.ids)
  VariableFeatures(object) <- features.integration
  object <- integrateData(object, structured_mat, assay = assay, new.assay.name = new.assay.name, nDims.scMC = nDims.scMC, lambda = lambda)
  return(object)
}




