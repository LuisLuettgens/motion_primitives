#include "matrix.h"

// fall M=1
template<int N, int M>
double Matrix<N,M>::norm() const {
	double ret = 0;
	for (int i=0;i<N;i++) {
		ret += D[i][0]*D[i][0];
	}
	return sqrt(ret);
}

template<int N, int M>
Matrix<M,N> Matrix<N,M>::Transpose() const {
	Matrix<M,N> A;
	for (int i=0;i<M;i++) {
		for (int j=0;j<N;j++) {
			A(i,j) = D[j][i];
		}
	}
	return A;
}





template<int N, int M>
void qr(Matrix<M,N> &A, Matrix<M,M> &Q) {

	for (int j=0;j<N;j++) {
		//		cout << "#########  j = " << j << endl;

		// x bestimmen
		Matrix<M,1> x;
		for (int i=0;i<M;i++) {
			x(i) = i<j?0:A(i,j);
		}
		//		cout << "x" << x;

		// Q bestimmen
		double xnorm = x.norm();
		double k = xnorm * (x(j)>0 ? -1. : 1.);
		double beta = 1./(xnorm*(abs(x(j))+xnorm));
		Matrix<M,1> u;
		for (int i=0;i<M;i++) {
			u(i) = i<j?0:x(i);
		}
		u(j) -= k;
		Matrix<M,M> Qi;

		//		cout << "J" << j << " " << beta << endl;
		//		cout << ":" << u << endl;
		for (int i=0;i<M;i++) {
			for (int k=0;k<M;k++) {
				Qi(i,k) = ((i==k)?1:0) - beta*u(i)*u(k);
			}
		}

		//		cout << "Q" << Q;
		//		cout << "qi" << Qi;
		A = Qi*A;
		//		cout << "Ai" << A;
		if (j==0)
			Q = Qi;
		else
			Q= Q*Qi;


		//		cout << "Q" << Q;
	}
}


template<int N, int M>
Matrix<M,1> Matrix<N,M>::Rueck(const Matrix<N,1> &c) const {
	Matrix<M,1> x;

	// Rx = c
	for (int i=N-1;i>=0;i--) {
		x(i) = c(i);
		for (int j=i+1;j<M;j++) {
			x(i) -= D[i][j] * x(j);
		}
		x(i) /= D[i][i];
	}
	return x;
}

template <int N, int M>
ostream &operator<<(ostream &os, const Matrix<N,M> &m) {
	os.setf(ios::fixed)
		;
	os << setprecision(4);
	for (int i=0;i<N;i++) {
		for (int j=0;j<M;j++) {
			os << setw(14) << m(i,j);
		}
		os << endl;
	}
	return os;
}

template <int N, int M, int L>
Matrix<N,L> operator*(const Matrix<N,M> &a, const Matrix<M,L> &b) {

	Matrix<N,L> A;
	for (int i=0;i<N;i++) {
		for (int j=0;j<L;j++) {
			A(i,j)=0;
			for (int k=0;k<M;k++) {
				A(i,j) += a(i,k) * b(k,j);
			}
		}
	}
	return A;

}

/*
void test() {

	Matrix<4,1> a;
	Matrix<4,4> c,b;
	b.Transpose();
	a.norm();
	qr(b,c);
	c.Rueck(a);

	std::cout << a<<c;
	b*c;
	c*a;
}
*/