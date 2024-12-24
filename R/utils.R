#' A function that performs basic utilities about the data based on lolR (Eric).

spp.utils.info <- function(X, Y, robust=FALSE, ...) {
  # Load necessary library
  require(robustbase)

  # Identify unique class labels
  unique_labels <- as.vector(sort(unique(Y)))

  num_samples <- nrow(X)
  num_features <- ncol(X)
  # Number of classes
  num_classes <- length(unique_labels)
  # compute the fraction of each class for prior p
  priors = sapply(unique_labels, function(label) sum(Y == label)/num_samples)
  if (!robust) {
    centroids <- as.matrix(array(sapply(unique_labels, function(label) colMeans(X[Y==label,,drop=FALSE])), dim=c(num_features, num_classes)))
  } else {
    # robust estimator of the mean is the median
    centroids <- as.matrix(array(sapply(unique_labels, function(label) colMedians(X[Y==label,,drop=FALSE])), dim=c(num_features, num_classes)))
  }
  return(list(n=num_samples, d=num_features, ylabs=unique_labels, priors=priors, centroids=centroids, K=num_classes))
}



spp.utils.deltas <- function(centroids, priors, ...) {
  d <- dim(centroids)[1]; K <- length(priors)
  # compute the rank-K difference space as deltas(i) = mus(i) - mus(0) where the mus are ordered
  # by decreasing prior
  deltas <- array(0, dim=c(d, K))
  srt_prior <- sort(priors, decreasing=TRUE, index.return=TRUE)$ix
  deltas[,1] <- centroids[,srt_prior[1]]
  for (i in 2:K) {
    deltas[,i] <- deltas[,1] - centroids[,srt_prior[i]]
  }
  return(deltas[,2:K,drop=FALSE])
}


spp.utils.ohe <- function(Y) {
  # Get the unique class labels and number of classes
  class_labels <- unique(Y)
  num_classes <- length(class_labels)
  num_samples <- length(Y)

  # Initialize the one-hot encoding matrix with zeros
  one_hot_encoded <- matrix(0, nrow = num_samples, ncol = num_classes)

  # Assign 1 where the sample belongs to a particular class
  for (i in seq_along(class_labels)) {
    one_hot_encoded[Y == class_labels[i], i] <- 1
  }

  # Assign column names to the one-hot encoded matrix
  colnames(one_hot_encoded) <- class_labels

  return(one_hot_encoded)
}
