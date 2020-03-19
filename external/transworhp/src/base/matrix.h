#ifndef matrix_h
#define matrix_h
#include "defines.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

template <int N, int M>
class Matrix {
public:
	double &operator()(int x, int y=0) {
		return D[x][y];
	}
	double operator()(int x, int y=0) const {
		return D[x][y];
	}

	double norm() const; // falls M=1
	Matrix<M,N> Transpose() const;
	Matrix<M,1> Rueck(const Matrix<N,1> &c) const;

private:
	double D[N][M];
};

template<int N, int M>
void qr(Matrix<M,N> &A, Matrix<M,M> &Q);

template <int N, int M>
ostream &operator<<(ostream &os, const Matrix<N,M> &m);

template <int N, int M, int L>
Matrix<N,L> operator*(const Matrix<N,M> &a, const Matrix<M,L> &b);

template class Matrix<4,1>
;
template class Matrix<4,4>
;


#endif
