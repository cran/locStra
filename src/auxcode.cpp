#include <R.h>
#include <Rcpp.h>
using namespace Rcpp;
#include <RcppEigen.h>
using namespace Eigen;

/*
using Eigen::Map;
using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::MatrixXi;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
*/



// ************************************************************************************
// standard matrix methods
// ************************************************************************************

// power method to compute the largest eigenvalue, set v=0 for automatic rowStarting value
// [[Rcpp::export]]
VectorXd powerMethodCpp(MatrixXd X, VectorXd v, double eps=1e-6, int maxiter=100) {
	if(X.rows()!=X.cols()) {
		Rcpp::stop("powerMethod requires square numeric matrix.");
	}
	
	if(v.size()<X.rows()) {
		v = VectorXd::Zero(X.rows());
		v = v.array()+1;
	}
	
	if(v.size()!=X.rows()) {
		Rcpp::stop("powerMethod requires X and v to be compatible.");
	}
	
	if(eps<=0) {
		eps = 1e-6;
	}
	
	VectorXd v_old(v);
	VectorXd v_new;
	VectorXd absdiff;
	for(int steps=0; steps<maxiter; steps++) {
		v_new = X * v_old;
		v_new.normalize();
		absdiff = v_new.array().abs() - v_old.array().abs();
		if(absdiff.norm() <= eps) {
			break;
		}
		v_old = v_new;
		steps = steps + 1;
	}
	return v_new;
}



// ************************************************************************************
// methods for dense matrices
// ************************************************************************************

// cpp function for computing the covariance matrix "cov" in R
// [[Rcpp::export]]
MatrixXd covCpp(MatrixXd X) {
	// scale: subtract column means from columns
	MatrixXd s = X.rowwise() - X.colwise().sum()*1/X.rows();
	// return scaled cross product
	return s.transpose() * s * 1/(X.rows()-1);
}



// fast vectorized computation of the Jaccard matrix
// [[Rcpp::export]]
MatrixXd jaccardMatrixCpp2(MatrixXd X) {
	VectorXd colsums = X.colwise().sum();
	MatrixXd shared = X.transpose() * X;
	MatrixXd total = ((shared.rowwise() - colsums.transpose()).colwise() - colsums).cwiseAbs();
	return shared.cwiseQuotient(total);
}



// get submatrix of dense matrix, specified by 0/1 selectors "cols" and "rows"
MatrixXd denseSubMatrix(MatrixXd X, VectorXd rows, VectorXd cols) {
	MatrixXd res = MatrixXd::Zero(int(rows.sum()), int(cols.sum()));
	int rowcounter = 0;
	int colcounter;
	for(int i=0; i<X.rows(); i++) {
		if(int(rows(i))==1) {
			colcounter = 0;
			for(int j=0; j<X.cols(); j++) {
				if(int(cols(j))==1) {
					res(rowcounter,colcounter) = X(i,j);
					colcounter++;
				}
			}
			rowcounter++;
		}
	}
	return res;
}

// get subset of rows from indicator vector "v" for dense matrix
MatrixXd denseSubsetRows(MatrixXd X, VectorXd v) {
	MatrixXd res = MatrixXd::Zero(int(v.sum()),X.cols());
	int rowcounter = 0;
	for(int i=0; i<v.size(); i++) {
		if(int(v(i))==1) {
			res.row(rowcounter) = X.row(i);
			rowcounter++;
		}
	}
	return res;
}

// cpp version of the s-matrix computation for dense input matrix
// [[Rcpp::export]]
MatrixXd calculateSMatrixDenseCpp(MatrixXd X, bool Djac=false, bool phased=false, int minVariants=5) {
	double numAlleles = X.cols();
	if(!phased) numAlleles *= 2.0;
	VectorXd sumVariants = X.rowwise().sum();
	
	for(int i=0; i<sumVariants.size(); i++) {
		if(sumVariants(i)>0.5*numAlleles) {
			if(phased) X.row(i) = 1 - X.row(i).array();
			else X.row(i) = 2 - X.row(i).array();
		}
	}
	
	sumVariants = X.rowwise().sum();
	for(int i=0; i<sumVariants.size(); i++) {
		if(sumVariants(i)>minVariants) sumVariants(i) = 1;
		else sumVariants(i) = 0;
	}
	MatrixXd Y = denseSubsetRows(X,sumVariants);
	
	VectorXd sumFilteredVariants = Y.rowwise().sum();
	double totalPossiblePairs = numAlleles*(numAlleles-1)/2.0;
	VectorXd totalPairs = sumFilteredVariants.array()*(sumFilteredVariants.array()-1)*0.5;
	VectorXd weights = totalPossiblePairs*totalPairs.array().inverse();
	
	MatrixXd s_matrix_numerator;
	if(!Djac) {
		MatrixXd Z = MatrixXd::Zero(Y.rows(),Y.cols());
		// multiply columns of Y with "weights" vector and store in Z
		for(int i=0; i<Y.cols(); i++) {
			Z.col(i) = Y.col(i).cwiseProduct(weights);
		}
		s_matrix_numerator = (Z.transpose() * Y).cast<double>();
	}
	else {
		s_matrix_numerator = (Y.transpose() * Y).cast<double>();
	}
	MatrixXd s_matrix_hap = s_matrix_numerator * (1.0/Y.rows());
	
	MatrixXd s_matrix_dip;
	if(phased) {
		// alternative binary vectors for rows and columns
		VectorXd rowsTF = VectorXd::Zero(s_matrix_hap.rows());
		for(int i=0; i<rowsTF.size(); i++) rowsTF(i) = (i+1)%2;
		VectorXd colsTF = VectorXd::Zero(s_matrix_hap.cols());
		for(int i=0; i<colsTF.size(); i++) colsTF(i) = (i+1)%2;
		s_matrix_dip = ( denseSubMatrix(s_matrix_hap,rowsTF,colsTF) + denseSubMatrix(s_matrix_hap,1-rowsTF.array(),colsTF) +
							denseSubMatrix(s_matrix_hap,rowsTF,1-colsTF.array()) + denseSubMatrix(s_matrix_hap,1-rowsTF.array(),1-colsTF.array()) )/4.0;
	} else {
		s_matrix_dip = s_matrix_hap/4.0;
	}
	return s_matrix_dip;
}



// genomic relationship matrix (classic and robust version of the grm) as defined in
// "Efficient Estimation of Realized Kinship from Single Nucleotide Polymorphism Genotypes" (Wang, Sverdlov, Thompson, 2017)
// [[Rcpp::export]]
MatrixXd grmDenseCpp3(MatrixXd X, bool robust=true) {
	// compute population frequencies across rows
	VectorXd p = X.rowwise().mean();
	VectorXd counterp = 1.0-p.array();
	VectorXd q = 4.0*p.cwiseProduct(counterp);
	// compute grm
	X = X.colwise() - 2*p;
	if(robust) return X.transpose() * X * 1.0/q.sum();
	else return X.transpose() * (X.array().colwise() * q.array().inverse()).matrix() * 1/double(X.rows());
}



// ************************************************************************************
// sparse matrix methods
// ************************************************************************************

// all inputs into sparse functions is a dense n*3 index matrix which contains non-zero entries as (x,y,value) tripels per row

// cpp function to write a sparse matrix
SparseMatrix<int> triplesToMatrix(MatrixXi T, int nrows, int ncols, int rowStart=0, int rowEnd=-1) {
	if(rowEnd==-1) rowEnd = nrows-1;
	SparseMatrix<int> X(nrows,ncols);
	for(int i=0; i<T.rows(); i++) {
		if( (T(i,0)>=rowStart) && (T(i,0)<=rowEnd) ) X.insert(T(i,0)-rowStart,T(i,1)) = T(i,2);
	}
	return X;
}

// copy sparse int matrix
SparseMatrix<int> copySparse(SparseMatrix<int> X) {
	SparseMatrix<int> Y(X.rows(),X.cols());
	for(int k=0; k<X.outerSize(); k++) {
		for(SparseMatrix<int>::InnerIterator it(X,k); it; ++it) {
			Y.insert(it.row(),it.col()) = it.value();
		}
	}
	return Y;
}

// copy subset of rows: vector v has to be 0/1 selector for rows
SparseMatrix<int> subsetRows(SparseMatrix<int> X, VectorXd v) {
	// assign new row numbers
	VectorXi w = VectorXi::Zero(v.size());
	int counter = 0;
	for(int i=0; i<v.size(); i++) {
		if(int(v(i))==1) {
			w(i) = counter;
			counter++;
		}
	}
	// now copy only rows in v
	SparseMatrix<int> Y(counter,X.cols());
	for(int k=0; k<X.outerSize(); k++) {
		for(SparseMatrix<int>::InnerIterator it(X,k); it; ++it) {
			if(int(v(it.row()))==1) Y.insert(w(it.row()),it.col()) = it.value();
		}
	}
	return Y;
}

VectorXd sparseRowSums(SparseMatrix<int> X) {
	VectorXd v = VectorXd::Zero(X.rows());
	for(int k=0; k<X.outerSize(); k++) {
		for(SparseMatrix<int>::InnerIterator it(X,k); it; ++it) {
			v[it.row()] += it.value();
		}
	}
	return v;
}

VectorXd sparseColSums(SparseMatrix<int> X) {
	VectorXd v = VectorXd::Zero(X.cols());
	for(int k=0; k<X.outerSize(); k++) {
		for(SparseMatrix<int>::InnerIterator it(X,k); it; ++it) {
			v[it.col()] += it.value();
		}
	}
	return v;
}



// cpp function for computing the covariance matrix "cov" in R, input is n*3 index matrix
// [[Rcpp::export]]
MatrixXd covCpp_sparse(MatrixXi T, int nrows, int ncols, int rowStart=0, int rowEnd=-1) {
	SparseMatrix<int> X = triplesToMatrix(T,nrows,ncols,rowStart,rowEnd);
	// return cov matrix
	VectorXd w = sparseColSums(X);
	VectorXd v = w * 1/X.rows();
	MatrixXd temp = (X.transpose()*X).cast<double>();
	return 1/(double(X.rows()-1)) * ( temp - w*v.transpose() - v*w.transpose() + X.rows()*v*v.transpose() );
}



// faster vectorized Jaccard matrix computation
// [[Rcpp::export]]
MatrixXd jaccardMatrixCpp4_sparse(MatrixXi T, int nrows, int ncols) {
	SparseMatrix<int> X = triplesToMatrix(T,nrows,ncols);
	VectorXd colsums = sparseColSums(X);
	MatrixXd shared = (X.transpose() * X).cast<double>();
	MatrixXd total = ((shared.rowwise() - colsums.transpose()).colwise() - colsums).cwiseAbs();
	return shared.cwiseQuotient(total);
}



// cpp code to calculate the S-matrix: function allows to perform computation directly on subpart [rowStart,rowEnd] of sparse matrix
// [[Rcpp::export]]
MatrixXd calculateSMatrixCpp(MatrixXi T, int nrows, int ncols, bool Djac=false, bool phased=false, int minVariants=5, int rowStart=0, int rowEnd=-1) {
	SparseMatrix<int> X = triplesToMatrix(T,nrows,ncols,rowStart,rowEnd);
	
	double numAlleles = X.cols();
	if(!phased) numAlleles *= 2.0;
	VectorXd sumVariants = sparseRowSums(X);
	
	for(int k=0; k<X.outerSize(); k++) {
		for(SparseMatrix<int>::InnerIterator it(X,k); it; ++it) {
			if(sumVariants(it.row())>0.5*numAlleles) {
				if(phased) X.coeffRef(it.row(),it.col()) = 1 - it.value();
				else X.coeffRef(it.row(),it.col()) = 2 - it.value();
			}
		}
	}
	
	sumVariants = sparseRowSums(X);
	for(int i=0; i<sumVariants.size(); i++) {
		if(sumVariants(i)>minVariants) sumVariants(i) = 1;
		else sumVariants(i) = 0;
	}
	SparseMatrix<int> Y = subsetRows(X,sumVariants);
	
	VectorXd sumFilteredVariants = sparseRowSums(Y);
	double totalPossiblePairs = numAlleles*(numAlleles-1)/2.0;
	VectorXd totalPairs = sumFilteredVariants.array()*(sumFilteredVariants.array()-1)*0.5;
	VectorXd weights = totalPossiblePairs*totalPairs.array().inverse();
	
	MatrixXd s_matrix_numerator;
	if(!Djac) {
		SparseMatrix<int> Z = copySparse(Y);
		// multiply columns of Z with "weights" vector
		for(int k=0; k<Z.outerSize(); k++) {
			for(SparseMatrix<int>::InnerIterator it(Z,k); it; ++it) {
				Z.coeffRef(it.row(),it.col()) = it.value()*weights(it.row());
			}
		}
		s_matrix_numerator = (Z.transpose() * Y).cast<double>();
	}
	else {
		s_matrix_numerator = (Y.transpose() * Y).cast<double>();
	}
	MatrixXd s_matrix_hap = s_matrix_numerator * (1.0/Y.rows());
	
	MatrixXd s_matrix_dip;
	if(phased) {
		// alternative binary vectors for rows and columns
		VectorXd rowsTF = VectorXd::Zero(s_matrix_hap.rows());
		for(int i=0; i<rowsTF.size(); i++) rowsTF(i) = (i+1)%2;
		VectorXd colsTF = VectorXd::Zero(s_matrix_hap.cols());
		for(int i=0; i<colsTF.size(); i++) colsTF(i) = (i+1)%2;
		s_matrix_dip = ( denseSubMatrix(s_matrix_hap,rowsTF,colsTF) + denseSubMatrix(s_matrix_hap,1-rowsTF.array(),colsTF) +
							denseSubMatrix(s_matrix_hap,rowsTF,1-colsTF.array()) + denseSubMatrix(s_matrix_hap,1-rowsTF.array(),1-colsTF.array()) )/4.0;
	} else {
		s_matrix_dip = s_matrix_hap/4.0;
	}
	return s_matrix_dip;
}



// sparse version of the cpp code to compute the (classic and robust) grm matrix
// [[Rcpp::export]]
MatrixXd grmCpp_sparse(MatrixXi T, int nrows, int ncols, bool robust=true) {
	SparseMatrix<double> X(nrows,ncols);
	for(int i=0; i<T.rows(); i++) {
		X.insert(T(i,0),T(i,1)) = T(i,2);
	}
	// compute population frequencies across rows
	VectorXd p = VectorXd::Zero(X.rows());
	for(int k=0; k<X.outerSize(); k++) {
		for(SparseMatrix<double>::InnerIterator it(X,k); it; ++it) {
			p(it.row()) += it.value();
		}
	}
	p = p * 1/double(X.cols());
	VectorXd counterp = 1.0-p.array();
	VectorXd q = 4.0*p.cwiseProduct(counterp);
	// compute grm
	VectorXd twop = 2*p;
	if(robust) {
		MatrixXd temp = (X.transpose()*X).cast<double>();
		return ( ( (temp.colwise() - X.transpose()*twop).rowwise() - twop.transpose()*X ).array() + twop.dot(twop) ) * 1.0/q.sum();
	}
	else {
		VectorXd invq = q.array().inverse();
		SparseMatrix<double> Y(nrows,ncols);
		for(int i=0; i<T.rows(); i++) {
			Y.insert(T(i,0),T(i,1)) = T(i,2) * invq(T(i,0));
		}
		MatrixXd temp = (X.transpose() * Y).cast<double>();
		return ( ( (temp.colwise() - X.transpose()*(twop.cwiseProduct(invq))).rowwise() - twop.transpose()*Y ).array() + twop.dot(twop.cwiseProduct(invq)) ) * 1.0/double(X.rows());
	}
}
