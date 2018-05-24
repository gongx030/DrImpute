#' @importFrom stats prcomp cor dist median kmeans
#' @importFrom Rcpp evalCpp
#' @importFrom parallel mclapply
#' @importFrom irlba prcomp_irlba
#' @importFrom FNN knnx.index
#' @import Matrix 
#' @useDynLib DrImpute
NULL

#' DrImpute
#'
#' Imputing dropout events in single-cell RNA-sequencing data.
#'
#' @param X Gene expression matrix (gene by cell). 
#' @param ks Number of cell clustering groups. Default set to ks = 10:15. 
#' @param dists Distribution matrices to use. Default is set to c("spearman", "pearson"). "eucleadian" can be added as well. 
#' @param fast Whether or not using approximation for large scale scRNA-seq data (default: FALSE)
#' @param dropout.probability.threshold This argument determine the dropout events for imputation.  If this argument is zero, all
#'				zero events will be considered for imputation; otherwise, only the zero entries that have greater probability than the specified threshold
#'				will be imputed (default: 0)
#' @param n.dropout Number of dropout entries for evaluting the dropout probability (default: 10000)
#' @param n.background Number of background entries for evaluting the dropout probability (default: 10000)
#' @param mc.cores Number of CPU cores for the fast method (default: 1)
#
#' @return A matrix object 
#' 
#' @export
#
#' @author Wuming Gong and Il-Youp Kwak
#'
#' @references
#' Wuming Gong, Il-Youp Kwak, Pruthvi Pota, Kaoko Koyano-Nakagawa and Daniel J. Garry (2017+)
#' DrImpute: Imputing dropout eveents in single cell RNA sequencing data
#'
#' @examples
#'
#' library(scDatasets)
#' library(SummarizedExperiment)
#' data(usoskin)
#' X <- assays(usoskin)$count
#' X <- preprocess(X, min.expressed.gene = 0)
#' X.log <- log(X + 1)
#' set.seed(1)
#' X.imp <- DrImpute(X.log)
#' 
DrImpute <- function(X, ks = 10:15, dists = c('spearman', 'pearson'), fast = FALSE, dropout.probability.threshold = 0, n.dropout = 10000, n.background = 10000, mc.cores = 1){

	N <- nrow(X)	# number of genes
	M <- ncol(X) 	# number of cells
	cat(sprintf('[%s] number of cells: %d\n', Sys.time(), M))
	cat(sprintf('[%s] number of genes: %d\n', Sys.time(), N))

	is.zero <- as(X == 0, 'lgeMatrix')	# dense logical matrix for zero and weakly detected events
	is.expressed <- as(X > log(1 + 5), 'lgeMatrix') # dense logical matrix for robustly expressed events
	cat(sprintf('[%s] %% zero events (X=0): %.1f%%\n', Sys.time(), sum(is.zero) / (N * M) * 100))
	cat(sprintf('[%s] %% robustly expressed events (X>5): %.1f%%\n', Sys.time(), sum(is.expressed) / (N * M) * 100))
	cat(sprintf('[%s] fast imputation: %s\n', Sys.time(), fast))

	if (fast){

		X <- as(X, 'dgCMatrix')
		if (dropout.probability.threshold > 0){
			dropout <- sample(Matrix::which(is.expressed), n.dropout)						# randomly sampled non-zero events
			x.dropout <- X[dropout]	# save the original expression values of artificial zeros
			X[dropout] <- 0	# generate artificial zeros
			background <- sample(Matrix::which(!is.zero), n.background)						# randomly sampled background zero or weekly detected events

			Xe <- estimate.expression(X, n = c(dropout, background), ks = ks, dists = dists, mc.cores = mc.cores)
			y <- rep(c(TRUE, FALSE), c(n.dropout, n.background))
			model <- glm(y ~ ., data.frame(y = y, Xe), family = 'binomial')
			X[dropout] <- x.dropout	# recover the original input

			Xe <- estimate.expression(X, n = is.zero, ks = ks, dists = dists, mc.cores = mc.cores)
			yp <- predict(model, data.frame(y = rep(FALSE, sum(is.zero)), Xe), type = 'response')	# predicting the dropout events vs. background zeros
			is.dropout <- yp > dropout.probability.threshold
			cat(sprintf('[%s] %% of predicted dropout among all zero events: %.1f%%\n', Sys.time(), sum(is.dropout) / length(is.dropout) * 100))
		}else{
			# imputing every zero event
			Xe <- estimate.expression(X, n = is.zero, ks = ks, dists = dists, mc.cores = mc.cores)
			is.dropout <- rep(TRUE, sum(is.zero))
		}

		xe <- rep(0, sum(is.zero))
		xe[is.dropout] <- Matrix::rowMeans(Xe[is.dropout, , drop = FALSE])
		X[is.zero] <- xe
	
	}else{
		cls <- getCls(X = X, ks = ks, dists = dists)
		czp <- sum(X == 0) / (N * M)
		Xe <- imp0clC(as.matrix(X), cls)
		rownames(Xe) <- rownames(X)
		colnames(Xe) <- colnames(X)
		X <- Xe
	}


	is.zero <- as(X == 0, 'lgeMatrix')	# dense logical matrix for zero indicators
	cat(sprintf('[%s] %% zero events after imputation: %.1f%%\n', Sys.time(), sum(is.zero) / (N * M) * 100))
	X

} # end of DrImpute


#' getCls
#'
#' get base clustering results using SC3 based clustering methods.
#' Similarity matrix constructed using "pearson", "spearman" or "euclidean". K-means clustering is performed on first few number of principal components of similarity matrix.
#' 
#' @param X Log transformed gene expression matrix (Gene by Cell). 
#' @param ks Number of cell clustering groups. Default set to ks = 10:15. 
#' @param dists Distribution matrices to use. Default is set to c("spearman", "pearson"). "euclidean" can be added as well. 
#' @param dim.reduc.prop Proportion of principal components to use for K-means clustering.
#'
#' @return A matrix object, Each row represent different clustering results.
#'
#' @author Il-Youp Kwak
#'
#' @references
#' Wuming Gong, Il-Youp Kwak, Pruthvi Pota, Kaoko Koyano-Nakagawa and Daniel J. Garry (2017+)
#' DrImpute: Imputing dropout eveents in single cell RNA sequencing data
getCls <- function(X, ks = 10:15, dists = c('spearman', 'pearson'), dim.reduc.prop = 0.05){

	N = dim(X)[2]
	ddim = round(N * dim.reduc.prop)

	Ds <- NULL
	i = 1
	for( d in dists) {
		if(d == "spearman") {
			Ds[[i]] <- as.matrix(1 - cor(X, method = "spearman"))
			i = i + 1;
		} else if (d == "pearson") {
			Ds[[i]] <- as.matrix(1 - cor(X, method = "pearson"))
			i = i + 1;
		} else if (d == "euclidean") {
			Ds[[i]] <- as.matrix(dist(t(X), method = "euclidean"))
			i = i + 1;
		}
	}

	cls <- NULL
	for(k in ks) {
		for(d in 1:length(Ds)) {
			cls <- rbind(cls, kmeans(prcomp(Ds[[d]], center = TRUE, scale. = TRUE)$rotation[,1:ddim], centers = k, iter.max = 1e+09, nstart = 1000)$cluster )
		}
	}
	cls

} # end of getCls


#' estimate.expression
#'
#' Estimate the expected expression levels for entries in the gene by cell expression matrix 
#'
#' @param X Log transformed gene expression matrix (Gene by Cell). 
#' @param n The entries for imputation. (default: NULL)
#' @param ks Number of cell clustering groups. (default: 10:15)
#' @param dists Distribution matrices to use. (default: c("spearman", "pearson"))
#' @param dim.reduc.prop Proportion of principal components to use for K-means clustering. (default: 0.05)
#' @param max.dim Maximum dimensions for PCA (default: 100)
#' @param pc.cdr.cc.cutoff The PC's which correlation coefficient with CDR (cellular detection rate) will be removed. (default: 1)
#' @param mc.cores Number of CPU cores (default: 1)
#'
#' @return A matrix object
#'
#' @author Wuming Gong
#'
#' @references
#' Wuming Gong, Il-Youp Kwak, Pruthvi Pota, Kaoko Koyano-Nakagawa and Daniel J. Garry (2017+)
#' DrImpute: Imputing dropout eveents in single cell RNA sequencing data
estimate.expression <- function(X, n = NULL, ks = 10:15, dists = c('spearman', 'pearson'), dim.reduc.prop = 0.05, max.dim = 100, pc.cdr.cc.cutoff = 1, mc.cores = 1){

	N <- nrow(X)	# number of genes
	M <- ncol(X) 	# number of cells

	ddim <- round(M * dim.reduc.prop)	# number of PCs for K-means clustering
	cdr <- Matrix::colSums(X > 0)	# cell-wise cell detection rate (CDR)
	is.zero <- as(X == 0, 'lgeMatrix')  # dense logical matrix for zero indicators

	if (ddim > max.dim){
		ddim <- max.dim 
	}

	pc.list <- lapply(dists, function(d){
		pc <- pca.dist.matrix(X, k = ddim, dist.method = d, mc.cores = mc.cores)	# PCA of the distance matrix 
		cc <- cor(pc, cdr)
		exclude <- abs(cc) > pc.cdr.cc.cutoff
		pc <- pc[, !exclude, drop = FALSE]	# remove the PC's that are highly correlated with the CDR
		pc
	})

	# clustering single cells and estimate expression based on clustering results
	param <- expand.grid(d = 1:length(pc.list), k = ks)
	Xp <- do.call('cbind', mclapply(1:nrow(param), function(i){
		d <- param[i, 'd']	# distance type
		k <- param[i, 'k']	# number of clusters
		mem <- kmeans2(pc.list[[d]], centers = k)$cluster	# k-means clustering
		M2C <- sparseMatrix(i = 1:M, j = mem, dims = c(M, k))	# cell ~ cluster member
		Xe <- (X %*% M2C) / ((!is.zero) %*% M2C)  # average expression levels of non-zero entries in each cluster
	  Xe[is.na(Xe) | is.infinite(Xe)] <- 0  # the genes that are all zeros in each cluster
	  Xe <- Xe %*% t(M2C) # expected expression levels, gene ~ cell
	  Xe[n]
	}, mc.cores = mc.cores))
	Xp

} # end of estimate.expression


#' pca.dist.matrix
#'
#' PCA of distance matrix
#'
#' @param X Log transformed gene expression matrix (Gene by Cell). 
#' @param k Number of PCs
#' @param dist.method Distance method (default: spearman)
#' @param mc.cores Number of CPU cores (default: 1)
#' @param batch.size Batch size of sampling based MDS approximation
#'
#' @return A matrix object
#'
#' @author Wuming Gong
#'
#' @references
#' Wuming Gong, Il-Youp Kwak, Pruthvi Pota, Kaoko Koyano-Nakagawa and Daniel J. Garry (2017+)
#' DrImpute: Imputing dropout eveents in single cell RNA sequencing data
pca.dist.matrix <- function(X, k, dist.method = 'spearman', mc.cores = 1, batch.size = 1000){

	M <- ncol(X)	# number of cells
	n.batch <- ceiling(M / batch.size)	# estimate how many batches
	groups <- sample(1:n.batch, M, replace = TRUE)	# randomly assign data into batches
	size <- table(factor(groups, 1:n.batch))	# number of samples per batch

	if (n.batch > 1)
		cat(sprintf('[%s] splitting data into %d batches\n', Sys.time(), n.batch))

	V.list <- mclapply(1:n.batch, function(b){
		D <- dist2(X[, groups == b], method = dist.method)
		t(prcomp2(D, n = k, center = TRUE, scale. = TRUE)$rotation)
	}, mc.cores = min(mc.cores, n.batch))	# PCA of cells within each batch

	if (n.batch > 1){
		m <- lapply(1:n.batch, function(b) sample(1:size[b], max(2, min(size[b], ceiling(batch.size / n.batch)))))	# for sampling a subset of cells from each batch
		m.align <- lapply(1:n.batch, function(b) which(groups == b)[m[[b]]])	# convert the local index to global index
		D <- dist2(X[, unlist(m.align)], method = dist.method)
		V.align <- t(prcomp2(D, n = k, center = TRUE, scale. = TRUE)$rotation)
		V.align <- lapply(split(1:ncol(V.align), list(rep(1:n.batch, sapply(m, length)))), function(i) V.align[, i, drop = FALSE])
		V.list <- mclapply(1:n.batch, function(b){	# map the local MDS to global MDS
			s <- svd(V.list[[b]][, m[[b]]])
			V.align[[b]] %*% s$v %*% diag(1 / s$d) %*% t(s$u) %*% V.list[[b]]
		}, mc.cores = mc.cores)
	}
	V <- matrix(0, k, M)
	for (b in 1:n.batch)
		V[, groups == b] <- V.list[[b]]
	t(V)
} # end of pca.dist.matrix


#' dist2
#'
#' Computing the distance matrix
#'
#' @param X Log transformed gene expression matrix (Gene by Cell). 
#' @param Distance method
#'
#' @return A matrix object
#'
#' @author Wuming Gong
#'
#' @references
#' Wuming Gong, Il-Youp Kwak, Pruthvi Pota, Kaoko Koyano-Nakagawa and Daniel J. Garry (2017+)
#' DrImpute: Imputing dropout eveents in single cell RNA sequencing data
dist2 <- function(X, method){
	if (method == 'spearman')
		D <- as.matrix(1 - cor(as.matrix(X), method = 'spearman'))
	else if (method == "pearson")
		D <- as.matrix(1 - cor(as.matrix(X), method = 'pearson'))
	else if (method == "euclidean") 
		D <- as.matrix(dist(t(as.matrix(X)), method = 'euclidean'))
	else
		stop(sprintf('unknown method: %s', method))
	D
} # end of dist2


#' prcomp2
#'
#' A wrapper for PCA
#'
#' @param X Log transformed gene expression matrix (Gene by Cell). 
#' @param n Number of PCs
#'
#' @return 
#'
#' @author Wuming Gong
#'
#' @references
#' Wuming Gong, Il-Youp Kwak, Pruthvi Pota, Kaoko Koyano-Nakagawa and Daniel J. Garry (2017+)
#' DrImpute: Imputing dropout eveents in single cell RNA sequencing data
prcomp2 <- function(X, n = NA, ...){
	if (is.na(n) || n > min(dim(X)) * 0.3){
		pr <- prcomp(X, ...)
		pr$sdev <- pr$sdev[1:n]
		pr$rotation <- pr$rotation[, 1:n, drop = FALSE]
		pr
	}else
		prcomp_irlba(X, n = n, ...)
} # end of prcomp2


#' kmeans2
#'
#' Mini-batch k-means
#'
#' @param x An input matrix
#' @param centers Number of centers
#' @param batch.size Batch size
#' @param iter.max Maximum iterations
#'
#' @references
#' Il-Youp Kwak, Wuming Gong, Kaoko Koyano-Nakagawa and Daniel J. Garry (2017+)
#' DrImpute: Imputing dropout eveents in single cell RNA sequencing data
#' 
#' D. Schully (2010) 
#' Web-scale k-means clustering
#'
kmeans2 <- function(x, centers, batch.size = 1000, iter.max = 1e+09){

	M <- nrow(x) # number of samples
	if (M < batch.size) 
		kmeans(x, centers = centers, iter.max = iter.max, nstart = 1000)
	else{
		# initialize each center with kmeans of randomly picked samples
		C <- kmeans(x[sample(1:M, batch.size), ], centers = centers, iter.max = iter.max, nstart = 1000)$centers
		v <- rep(1, centers)	# per-center counts
		
		iter <- 1
		while (iter <= iter.max){
			Cp <- C
			m <- sample(1:M, batch.size)	# samples randomly picked from x
			mem <- c(knnx.index(C, x[m, ], k = 1))	# center nearest to each sample
			S <- sparseMatrix(i = mem, j = 1:batch.size, dims = c(centers, batch.size))
			vc <- Matrix::rowSums(S)
			v <- v + vc	# upate per-center count
			w <- 1 / vc
			w[is.infinite(w)] <- 0
			S <- Diagonal(x = w) %*% S
			eta <- 1 / v	# per-center learning rate
			C <- Diagonal(x = 1 - eta) %*% C + Diagonal(x = eta) %*% S %*% x[m, ]	# take gradient step
			if (iter %% 10 == 0)
				if (norm(as.matrix(C - Cp), '2') < 1e-3)
					break
			iter <- iter + 1
		}
		list(cluster = c(knnx.index(C, x, k = 1)))
	}
} # end of kmeans2


#' softmax
#'
#' Computing y[i, j] = exp(x[i, j]) / sum(exp(x[, j]))
#'
#' @param x An input matrix
softmax <- function(x){
	x.max <- apply(x, 2, max)
	y <- log(Matrix::rowSums(exp(t(x) - x.max)) + .Machine$double.eps) + x.max	# log(sum(exp(x[, j])))
	y <- t(exp(t(x) - y))
} # end of softmax 

