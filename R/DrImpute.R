#' @importFrom stats prcomp cor dist median kmeans
#' @importFrom Rcpp evalCpp
#' @importFrom parallel mclapply
#' @importFrom irlba prcomp_irlba
#' @importFrom FNN knnx.index
#' @import Matrix 
#' @useDynLib DrImpute
NULL


#' A function for preprocessing gene expression matrix.
#'
#' Preprocess gene expression data 
#' 
#' @param X Gene expression matrix (Gene by Cell). 
#' @param min.expressed.gene Cell level filtering criteria. For a given cell, if the number of expressed genes are less than min.expressed.gene, we filter it out.  
#' @param min.expressed.cell Gene level filtering criteria. For a given gene, if the number of expressed cells are less than min.expressed.cell, we filter it out.  
#' @param max.expressed.ratio Gene level filtering criteria. For a given gene, if the ratio of expressed cells are larger than max.expressed.ratio, we filter it out.
#' @param normalize.by.size.effect Normaize using size factor.
#' @return Filtered gene expression matrix
#' @author Wuming Gong
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
#'
#' @seealso \code{\link{DrImpute}}
preprocess <- function(x, min.expressed.gene = 0, min.expressed.cell = 2, max.expressed.ratio = 1, normalize.by.size.effect = FALSE){

	if (class(x) == 'SummarizedExperiment')
		X <- assays(x)$count
	else if (class(x) == 'matrix')
		X <- x
	else if (is(x, 'sparseMatrix'))
		X <- x
	else
		stop(sprintf('unknown class(x): %s', class(x)))

	M <- ncol(X)
	N <- nrow(X)
	m <- Matrix::colSums(X > 1) >= min.expressed.gene	# cells that have at least min.expressed.gene expreseed genes
	n <- Matrix::rowSums(X > 1) <= max.expressed.ratio * M & Matrix::rowSums(X > 1) >= min.expressed.cell	# genes that are detected in at least min.expressed.cell or at most max.expressed.ratio cells
	if (normalize.by.size.effect){
	  sf <- apply((X[n, m] + 1) / exp(Matrix::rowMeans(log(X[n, m] + 1))), 2, median)
		X <- t(t(X[n, m]) / sf)
	}else
		X <- X[n, m]

	if (class(x) == 'SummarizedExperiment'){
		x <- x[n, m]
		assays(x)$count <- X
	}else if (class(x) == 'matrix'){
		x <- as.matrix(X)
	}else if (is(x, 'sparseMatrix')){
		x <- X
	}
	
	x
} # end of preprocess


#' DrImpute
#'
#' Imputing dropout events in single-cell RNA-sequencing data.
#'
#' @param X Gene expression matrix (gene by cell). 
#' @param ks Number of cell clustering groups. Default set to ks = 10:15. 
#' @param dists Distribution matrices to use. Default is set to c("spearman", "pearson"). "eucleadian" can be added as well. 
#' @param mc.cores Number of CPU cores (default: 1)
#' @param batch.size Batch size of sampling based MDS approximation; if batch.size is NA, the standard MDS is used. 
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
DrImpute <- function(X, ks = 10:15, dists = c('spearman', 'pearson'), batch.size = NA, mc.cores = 1){

	N <- nrow(X)	# number of genes
	M <- ncol(X) 	# number of cells
	cat(sprintf('[%s] number of cells: %d\n', Sys.time(), M))
	cat(sprintf('[%s] number of genes: %d\n', Sys.time(), N))

	is.zero <- as(X == 0, 'lgeMatrix')	# dense logical matrix for zero and weakly detected events
	is.expressed <- as(X > log(1 + 5), 'lgeMatrix') # dense logical matrix for robustly expressed events
	cat(sprintf('[%s] %% zero events (X=0): %.1f%%\n', Sys.time(), sum(is.zero) / (N * M) * 100))
	cat(sprintf('[%s] %% robustly expressed events (X>5): %.1f%%\n', Sys.time(), sum(is.expressed) / (N * M) * 100))
	if (!is.na(batch.size)){
		cat(sprintf('[%s] batch size: %d cells\n', Sys.time(), batch.size))
	}
	cat(sprintf('[%s] # of CPU cores: %d\n', Sys.time(), mc.cores))

	X <- estimate.expression(X, ks = ks, dists = dists, mc.cores = mc.cores, batch.size = batch.size)
	is.zero <- as(X == 0, 'lgeMatrix')	# dense logical matrix for zero indicators
	cat(sprintf('[%s] %% zero events after imputation: %.1f%%\n', Sys.time(), sum(is.zero) / (N * M) * 100))
	X

} # end of DrImpute


#' estimate.expression
#'
#' Estimate the expected expression levels for entries in the gene by cell expression matrix 
#'
#' @param X Log transformed gene expression matrix (Gene by Cell). 
#' @param ks Number of cell clustering groups. (default: 10:15)
#' @param dists Distribution matrices to use. (default: c("spearman", "pearson"))
#' @param dim.reduc.prop Proportion of principal components to use for K-means clustering. (default: 0.05)
#' @param max.dim Maximum dimensions for PCA (default: 100)
#' @param pc.cdr.cc.cutoff The PC's which correlation coefficient with CDR (cellular detection rate) will be removed. (default: 1)
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
estimate.expression <- function(X, ks = 10:15, dists = c('spearman', 'pearson'), dim.reduc.prop = 0.05, max.dim = 100, pc.cdr.cc.cutoff = 1, mc.cores = 1, batch.size = NA){

	X <- as(X, 'dgCMatrix')
	N <- nrow(X)	# number of genes
	M <- ncol(X) 	# number of cells

	ddim <- round(M * dim.reduc.prop)	# number of PCs for K-means clustering
	cdr <- Matrix::colSums(X > 0)	# cell-wise cell detection rate (CDR)
	is.zero <- as(X == 0, 'lgeMatrix')  # dense logical matrix for zero indicators

	if (ddim > max.dim){
		ddim <- max.dim 
	}

	pc.list <- lapply(dists, function(d){
		pc <- pca.dist.matrix(X, k = ddim, dist.method = d, mc.cores = mc.cores, batch.size = batch.size)	# PCA of the distance matrix 
		cc <- cor(pc, cdr)
		exclude <- abs(cc) > pc.cdr.cc.cutoff
		pc <- pc[, !exclude, drop = FALSE]	# remove the PC's that are highly correlated with the CDR
		pc
	})

	# clustering single cells and estimate expression based on clustering results
	param <- expand.grid(d = 1:length(pc.list), k = ks)
	cls <- do.call('rbind', mclapply(1:nrow(param), function(i){
		d <- param[i, 'd']	# distance type
		k <- param[i, 'k']	# number of clusters
		mem <- kmeans2(pc.list[[d]], centers = k)$cluster	# k-means clustering
		mem
	}, mc.cores = mc.cores))

	cat(sprintf('[%s] imputing zeros\n', Sys.time()))
	Xe <- imp0clC(as.matrix(X), cls)
	rownames(Xe) <- rownames(X)
	colnames(Xe) <- colnames(X)
	Xe	

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
pca.dist.matrix <- function(X, k, dist.method = 'spearman', mc.cores = 1, batch.size = NA){

	M <- ncol(X)	# number of cells

	if (is.na(batch.size))
		batch.size <- M

	n.batch <- ceiling(M / batch.size)	# estimate how many batches

	if (n.batch == 1){

		D <- dist2(X, method = dist.method)
		V <- t(prcomp2(D, n = k, center = TRUE, scale. = TRUE)$rotation)

	}else{

		groups <- sample(1:n.batch, M, replace = TRUE)	# randomly assign data into batches
		size <- table(factor(groups, 1:n.batch))	# number of samples per batch

		V.list <- mclapply(1:n.batch, function(b){
			D <- dist2(X[, groups == b], method = dist.method)
			t(cmdscale(D, k))
		}, mc.cores = min(mc.cores, n.batch))	# PCA of cells within each batch

		m <- lapply(1:n.batch, function(b) sample(1:size[b], max(2, min(size[b], ceiling(batch.size / n.batch)))))	# for sampling a subset of cells from each batch
		m.align <- lapply(1:n.batch, function(b) which(groups == b)[m[[b]]])	# convert the local index to global index
		D <- dist2(X[, unlist(m.align)], method = dist.method)
		V.align <- t(cmdscale(D, k))
		V.align <- lapply(split(1:ncol(V.align), list(rep(1:n.batch, sapply(m, length)))), function(i) V.align[, i, drop = FALSE])
		V.list <- mclapply(1:n.batch, function(b){	# map the local MDS to global MDS
			s <- svd(V.list[[b]][, m[[b]]])
			V.align[[b]] %*% s$v %*% diag(1 / s$d) %*% t(s$u) %*% V.list[[b]]
		}, mc.cores = mc.cores)

		V <- matrix(0, k, M)
		for (b in 1:n.batch)
			V[, groups == b] <- V.list[[b]]
	}
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
kmeans2 <- function(x, centers, batch.size = 10000, iter.max = 1e+09){

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



