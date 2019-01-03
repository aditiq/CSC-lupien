library(flexclust)

# Custom version of flexClust's kcca() function for k-centroid clustering.
kcca_Wouter <- function (x, k, family = kccaFamily("ejaccard"), control = NULL, stop_percentage=0.01) {
  MYCALL <- match.call()
  control <- as(control, "flexclustControl")
  N <- nrow(x)
  
  cat("Initialize cluster centers...")
  centers <- x[sample(1:nrow(x), k),];
  cluster <- integer(nrow(x));
  cat("done.\n");
  
  for (iter in 1:control@iter.max) {
    cat(paste("Iteration ", iter, ":\n", sep=""));
    clustold <- cluster
    
    cat("- Calculating distances to cluster centers...");
    distmat <- family@dist(x, centers)
    cat("done.\n");
    
    cat("- Obtain new cluster assignments...");
    cluster <- family@cluster(x, distmat = distmat)
    cat("done.\n");
    
    cat("- Update cluster centers...");
    #centers <- family@allcent(x, cluster = cluster, k = k)
    centers <- matrix(NA, nrow = k, ncol = ncol(x))
    for (n in unique(cluster)) {
      centers[n, ] <- colMeans(x[cluster == n, , drop=FALSE]);
    }
    
    centers <- centers[complete.cases(centers), , drop = FALSE]
    k <- nrow(centers)
    cat("done.\n");
    
    changes <- sum(cluster != clustold)
    cat(paste("--> CHANGES: ", changes, " (", round(changes/nrow(x)*100, 2), "%)\n\n", sep=""));
    
    # If less than 0.01% change, stop iterating
    if ((changes/nrow(x)) < (stop_percentage / 100)) { 
      break
    }
  }
  
  centers <- centers[complete.cases(centers), , drop = FALSE]
  distmat <- family@dist(x, centers)
  cluster <- family@cluster(x, distmat = distmat)
  
  z <- list(cluster = cluster, centers = centers, iter = iter, 
            converged = (iter < control@iter.max), 
            call = MYCALL, control = control, distance=distmat)
  
  z
}
