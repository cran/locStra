# **************************************************************************************************************************************
# Cpp implementation of local population stratification methods. This file provides Rcpp wrapper functions to compute:
# 
# 1) the windows to be used in function 'fullscan'.
# 2) the covariance matrix
# 3) the Jaccard similarity matrix
# 4) the sMatrix as in the 'Stego' package on https://github.com/dschlauch/stego
# 5) the classic and robust version of the genomic relationship matrix (grm)
# 6) the largest eigenvector via the power method (von Mises iteration algorithm), see https://en.wikipedia.org/wiki/Power_iteration
# 
# The main function is 'fullscan'. It can be used to compute global and local correlations in population stratification data.
# 
# The functions directly work on sparse matrix objects of the 'Matrix' class and never unpack the sparse matrix during computation.
# If the input is in fact a dense matrix, a separate 'dense' option is provided in all functions.
# Important note: For all functions the input matrix is always assumed to be oriented to contain the data for one individual per column.
# **************************************************************************************************************************************
#' @useDynLib locStra



# not exported
# auxiliary method to convert a sparse matrix into a dense matrix containing rows (i,j,x) encoding entry x at position (i,j) which are used as input to the C++ code
# input: sparse matrix 'm' of class 'Matrix'
# output: R matrix with rows (i,j,x) encoding entry x at position (i,j)
sparseToList <- function(m) {
	L <- summary(Matrix(m,sparse=TRUE))
	# cpp code counts from 0, usually summary() orders according to columns, thus transpose to obtain row sorting
	cbind(L$i-1,L$j-1,L$x)
}



#' Auxiliary function to generate a two-column matrix of window sizes to be used in the function 'fullscan'.
#' 
#' @param len The overall length of the data which is to be scanned in windows.
#' @param size The window size.
#' @param offset The offset of the generated windows (e.g., if \code{offset=1} then sliding window, if \code{offset=size} then blocks).
#' 
#' @return A two-column matrix of sliding windows, with one window per row defined through start and end value.
#' 
#' @importFrom Rdpack reprompt
#' 
#' @examples
#' library(locStra)
#' print(makeWindows(100,10,5))
#' 
#' @export
makeWindows <- function(len,size,offset) {
	if(len<size) return(matrix(c(1,len),nrow=1))
	w <- seq(1,len-(size-1),offset)
	wmatrix <- cbind(w,w+(size-1))
	lastw <- wmatrix[nrow(wmatrix),ncol(wmatrix)]
	if(lastw<len) {
		if(lastw+1<len) wmatrix <- rbind(wmatrix,c(lastw+1,len))
		else wmatrix[nrow(wmatrix),ncol(wmatrix)] <- wmatrix[nrow(wmatrix),ncol(wmatrix)]+1
	}
	unname(wmatrix)
}



#' Cpp implementation of the power method (von Mises iteration) to compute the largest eigenvector for a (sparse) input matrix.
#' 
#' @param m Symmetric matrix for which the largest eigenvector is sought.
#' @param initvector Optional vector compatible with the input matrix which serves as a starting value for the iteration. Default is zero.
#' 
#' @return The largest eigenvector of \code{m}.
#' 
#' @importFrom Rdpack reprompt
#' @references Richard von Mises and Hilda Pollaczek-Geiringer (1929). Praktische Verfahren der Gleichungsaufloesung. ZAMM Zeitschrift fuer Angewandte Mathematik und Mechanik, 9:152-164.
#' 
#' @examples
#' library(locStra)
#' m <- matrix(1:9,3)
#' print(powerMethod(m))
#' 
#' @export
powerMethod <- function(m,initvector=0) {
	powerMethodCpp(m,initvector)
}



#' Cpp implementation of a function to compute the covariance matrix for a (sparse) matrix. The function is equivalent to the R command 'cov' applied to matrices.
#' 
#' @param m A (sparse) matrix for which the covariance matrix is sought. The input matrix is assumed to be oriented to contain the data for one individual per column.
#' @param dense Flag to switch between purpose-built dense or sparse implementations. Default is \code{dense=FALSE}.
#' 
#' @return The covariance matrix of \code{m}.
#' 
#' @importFrom Rdpack reprompt
#' @references R Core Team (2014). R: A Language and Environment for Statistical Computing. R Foundation for Stat Comp, Vienna, Austria.
#' 
#' @examples
#' library(locStra)
#' library(Matrix)
#' m <- matrix(sample(0:1,15,replace=TRUE),ncol=3)
#' sparseM <- Matrix(m,sparse=TRUE)
#' print(covMatrix(sparseM))
#' 
#' @export
covMatrix <- function(m,dense=FALSE) {
	if(dense) {
		return(covCpp(as.matrix(m)))
	}
	else {
		return(covCpp_sparse(sparseToList(m),nrow(m),ncol(m)))
	}
}



#' Cpp implementation of the Jaccard similarity matrix computation for a (sparse) input matrix.
#' 
#' @param m A (sparse) matrix for which the Jaccard similarity matrix is sought. The input matrix is assumed to be oriented to contain the data for one individual per column.
#' @param dense Flag to switch between purpose-built dense or sparse implementations. Default is \code{dense=FALSE}.
#' 
#' @return The Jaccard matrix of \code{m}.
#' 
#' @importFrom Rdpack reprompt
#' @references Dmitry Prokopenko, Julian Hecker, Edwin Silverman, Marcello Pagano, Markus Noethen, Christian Dina, Christoph Lange and Heide Fier (2016). Utilizing the Jaccard index to reveal population stratification in sequencing data: a simulation study and an application to the 1000 Genomes Project. Bioinformatics, 32(9):1366-1372.
#' 
#' @examples
#' library(locStra)
#' library(Matrix)
#' m <- matrix(sample(0:1,15,replace=TRUE),ncol=3)
#' sparseM <- Matrix(m,sparse=TRUE)
#' print(jaccardMatrix(sparseM))
#' 
#' @export
jaccardMatrix <- function(m,dense=FALSE) {
	if(dense) {
		return(jaccardMatrixCpp2(as.matrix(m)))
	}
	else {
		return(jaccardMatrixCpp4_sparse(sparseToList(m),nrow(m),ncol(m)))
	}
}



#' Cpp implementation of the s-matrix function (which computes the weighted Jaccard similarity matrix) for a (sparse) input matrix as in the 'Stego' package on https://github.com/dschlauch/stego.
#' 
#' @param m A (sparse) matrix for which the s-matrix is sought. The input matrix is assumed to be oriented to contain the data for one individual per column.
#' @param dense Flag to switch between purpose-built dense or sparse implementations. Default is \code{dense=FALSE}.
#' @param phased Boolean flag to indicate if input matrix is phased. Default is \code{phased=FALSE}.
#' @param minVariants Integer cutoff value for minimal number of variants. Default is \code{minVariants=5}.
#' 
#' @return The s-matrix (the weighted Jaccard matrix) of \code{m}.
#' 
#' @importFrom Rdpack reprompt
#' @references Daniel Schlauch (2016). Implementation of the stego algorithm - Similarity Test for Estimating Genetic Outliers. https://github.com/dschlauch/stego
#' 
#' @examples
#' library(locStra)
#' library(Matrix)
#' m <- matrix(sample(0:1,15,replace=TRUE),ncol=3)
#' sparseM <- Matrix(m,sparse=TRUE)
#' print(sMatrix(sparseM))
#' 
#' @export
sMatrix <- function(m,dense=FALSE,phased=FALSE,minVariants=5) {
	Djac <- F
	if(dense) {
		return(calculateSMatrixDenseCpp(as.matrix(m),Djac,phased,minVariants))
	}
	else {
		return(calculateSMatrixCpp(sparseToList(m),nrow(m),ncol(m),Djac,phased,minVariants))
	}
}



#' Cpp implementation of the genomic relationship matrix (grm) for a (sparse) input matrix as defined in Yang et al. (2011).
#' 
#' @param m A (sparse) matrix for which the genomic relationship matrix is sought. The input matrix is assumed to be oriented to contain the data for one individual per column.
#' @param dense Flag to switch between purpose-built dense or sparse implementations. Default is \code{dense=FALSE}.
#' @param robust Flag to indicate if the classic (\code{robust=FALSE}) or robust (\code{robust=TRUE}) version of the grm matrix is sought. Default is \code{robust=TRUE}.
#' 
#' @return The genomic relationship matrix of \code{m}.
#' 
#' @importFrom Rdpack reprompt
#' @references Yang J, Lee SH, Goddard ME, Visscher PM (2011). GCTA: a tool for genome-wide complex trait analysis. Am J Hum Genet, 88(1):76-82.
#' 
#' @examples
#' library(locStra)
#' library(Matrix)
#' m <- matrix(sample(0:1,15,replace=TRUE),ncol=3)
#' sparseM <- Matrix(m,sparse=TRUE)
#' print(grMatrix(sparseM))
#' 
#' @export
grMatrix <- function(m,dense=FALSE,robust=TRUE) {
	if(dense) {
		return(grmDenseCpp3(as.matrix(m),robust))
	}
	else {
		return(grmCpp_sparse(sparseToList(m),nrow(m),ncol(m),robust))
	}
}



#' Main function: A full scan of the input data \code{m} using a collection of windows given by the two-column matrix \code{windows}. For each window, the data is processed using the function \code{matrixFunction} (this could be e.g. the \code{covMatrix} function), then the processed data is summarised using the function \code{summaryFunction} (e.g., the largest eigenvector computed with the function \code{powerMethod}), and finally the global and local summary scores (e.g., the largest eigenvectors) are compared using the function \code{comparisonFunction} (e.g., the vector correlation with R's function \code{cor}). The function returns a two-column matrix which contains per row the global (e.g., the correlation between global and local eigenvectors) and local (e.g., the correlation between the local eigenvector for the current window and the eigenvector for the last window) summary statistics for each window.
#' 
#' @param m A (sparse) matrix for which the full scan is sought. The input matrix is assumed to be oriented to contain the data for one individual per column.
#' @param windows A two-column matrix containing per column the windows on which the data is scanned. The windows can be overlapping. The windows can be computed using the function \code{makeWindows}.
#' @param matrixFunction Function on one matrix argument to process the data for each window (e.g., the covariance matrix).
#' @param summaryFunction Function on one argument to summarise the output of the function \code{matrixFunction} (e.g., the largest eigenvector).
#' @param comparisonFunction Function on two inputs to compute some kind of comparison measure for the output of the function \code{summaryFunction} (e.g., vector correlation, or matrix norm).
#' 
#' @return A two-column matrix containing per row the global and local summary statistics for each window.
#' 
#' @importFrom Rdpack reprompt
#' @references Dmitry Prokopenko, Julian Hecker, Edwin Silverman, Marcello Pagano, Markus Noethen, Christian Dina, Christoph Lange and Heide Fier (2016). Utilizing the Jaccard index to reveal population stratification in sequencing data: a simulation study and an application to the 1000 Genomes Project. Bioinformatics, 32(9):1366-1372.
#' 
#' @examples
#' library(locStra)
#' m <- matrix(sample(0:1,1000,replace=TRUE),ncol=10)
#' w <- makeWindows(nrow(m),10,10)
#' print(fullscan(m,w,covMatrix,powerMethod,cor))
#' 
#' @export
fullscan <- function(m,windows,matrixFunction,summaryFunction,comparisonFunction) {
	matrix_global <- matrixFunction(m)
	summary_global <- summaryFunction(matrix_global)
	last_summary <- rep(0,length(summary_global))
	
	# go through sliding window
	res <- matrix(0, nrow=nrow(windows), ncol=2)
	for(i in 1:nrow(windows)) {
		matrix_local <- matrixFunction(m[windows[i,1]:windows[i,2],])
		summary_local <- summaryFunction(matrix_local)
		res[i,] <- c(comparisonFunction(summary_global,summary_local), comparisonFunction(last_summary,summary_local))
	}
	return(res)
}
