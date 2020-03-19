#include "newdouble.h"
#include <iomanip>
#include "../base/vectortools.h"

using std::cout;
using std::endl;
using std::setw;
using std::vector;

namespace tw {

int MagicDouble::Derivate = 1;
int MagicDouble::Stats[10] = {0,0,0,0,0,0,0,0,0,0};

void MagicDouble::Info() {

	cout << " NEW      " << setw(10) << Stats[0] << endl;
	cout << " NO_DPNEW " << setw(10) << Stats[6] << endl;
	cout << " Const(d) " << setw(10) << Stats[1] << endl;
	cout << " Const()  " << setw(10) << Stats[2] << endl;
	
	cout << " dep      " << setw(10) << Stats[3] << endl;
	cout << " dep_also " << setw(10) << Stats[4] << endl;
	cout << " # dep    " << setw(10) << Stats[5] << endl;

	cout << " dep2     " << setw(10) << Stats[7] << endl;
	cout << " dep2also " << setw(10) << Stats[8] << endl;
	cout << " # dep2   " << setw(10) << Stats[9] << endl;

}

MagicDouble::MagicDouble() {

	Stats[2]++;
	val = 0;
}


MagicDouble::MagicDouble(double v) {

	Stats[1]++;
	val = v;

	MagicDeriv d;
	d.ref = this;
	d.diff = 1;
	depends.push_back(d);

}

MagicDouble::MagicDouble(MagicDouble &a, double v, double dv, double ddv) {

	Stats[1]++;
	val = v;

	{
		MagicDeriv d;
		d.ref = &a;
		d.diff = dv;
		depends.push_back(d);
	}

	{
		Magic2Deriv d;
		d.ref1 = &a;
		d.ref2 = &a;
		d.diff = ddv;
		depends2.push_back(d);
	}
}



void MagicDouble::link() {

	MagicDeriv d;
	d.ref = this;
	d.diff = 1;
	depends.push_back(d);

}


void MagicDouble::unlink() {
	
	depends.clear();
	depends2.clear();
}

void MagicDouble::COPY(MagicDouble *ret, const double *data, int n) {

	for (int i=0;i<n;i++) {
		ret[i].val = data[i];
	}

}

MagicDouble* MagicDouble::NEW(int n) {

	Stats[0]++;
	auto  ret = new MagicDouble[n];

	for (int i=0;i<n;i++) {

		MagicDeriv d;
		d.ref = &ret[i];
		d.diff = 1;
		ret[i].depends.push_back(d);

		Stats[1]++;
	}

	return ret;
}

MagicDouble* MagicDouble::NO_DEP_NEW(int n) {

	Stats[6]++;
	auto  ret = new MagicDouble[n];

	return ret;
}


void MagicDouble::DEBUG_(MagicDouble *ret, int n) {

	cout << "Block of size " << n << endl;
	for (int i=0;i<n;i++) {
		if (ret[i].depends.size()) {
			cout << "      " << ret[i] << endl;
		}
	}
	
}

const double* MagicDouble::DF(MagicDouble *d) const {

	auto fit = depends.begin();
	for (;fit!=depends.end();fit++) {
		if (fit->ref == d) {
			return &(fit->diff);
		}
	}

	return nullptr;

}

const double* MagicDouble::DDF(MagicDouble *d1, MagicDouble *d2) const {

	auto fit = depends2.begin();
	for (;fit!=depends2.end();fit++) {
		if (fit->ref1 == d1 && fit->ref2 == d2) {
			return &(fit->diff);
		}
	}

	return nullptr;

}


void MagicDouble::scale(double factor) {

	if (Derivate<1) return;

	auto it = depends.begin();
	for (;it!=depends.end();it++) {
		
		it->diff *= factor;
	}
}
void MagicDouble::scale2(double factor) {

	if (Derivate<2) return;

	auto it = depends2.begin();
	for (;it!=depends2.end();it++) {
		
		it->diff *= factor;
	}
}

void MagicDouble::depends_on(const vector<MagicDeriv> &other, double factor) {

	if (Derivate<1) return;

	Stats[3]++;
	auto it = other.begin();
	for (;it!=other.end();it++) {

		MagicDeriv d;
		d.ref = it->ref;
		d.diff = it->diff * factor;

		Stats[5]++;
		depends.push_back(d);
	}
}

void MagicDouble::depends_also_on(const vector<MagicDeriv> &other, double factor) {

	if (Derivate<1) return;

	Stats[4]++;
	auto it = other.begin();
	for (;it!=other.end();it++) {

		int found = 0;

		auto fit = depends.begin();
		for (;fit!=depends.end();fit++) {
			if (fit->ref == it->ref) {
				fit->diff += it->diff*factor;
				found = 1;
				break;
			}
		}

		if (found==0) {
			MagicDeriv d;
			d.ref = it->ref;
			d.diff = it->diff * factor;
			Stats[5]++;
			depends.push_back(d);
		}

	}
}




void MagicDouble::depends2_on(const vector<Magic2Deriv> &other, double factor) {

	if (Derivate<2) return;

	Stats[7]++;
	auto it = other.begin();
	for (;it!=other.end();it++) {

		Magic2Deriv d;
		d.ref1 = it->ref1;
		d.ref2 = it->ref2;

		d.diff = it->diff * factor;

		Stats[9]++;
		depends2.push_back(d);
	}
}

void MagicDouble::depends2_also_on(const vector<Magic2Deriv> &other, double factor) {

	if (Derivate<2) return;

	Stats[8]++;
	auto it = other.begin();
	for (;it!=other.end();it++) {

		int found = 0;

		auto fit = depends2.begin();
		for (;fit!=depends2.end();fit++) {
			
			if ((fit->ref1 == it->ref1) && (fit->ref2 == it->ref2)) {
				fit->diff += it->diff*factor;
				found = 1;
				break;
			}
		
		}

		if (found==0) {
			Magic2Deriv d;
			d.ref1 = it->ref1;
			d.ref2 = it->ref2;
	
			d.diff = it->diff * factor;
			Stats[9]++;
			depends2.push_back(d);
		}

	}
}




void MagicDouble::depends2_also_on(const vector<MagicDeriv> &other1, const vector<MagicDeriv> &other2, double factor) {

	if (Derivate<2) return;

	Stats[8]++;
	auto it1 = other1.begin();
	for (;it1!=other1.end();it1++) {

		auto it2 = other2.begin();
		for (;it2!=other2.end();it2++) {

			int found = 0;

			auto fit = depends2.begin();
			for (;fit!=depends2.end();fit++) {
						
				if ((fit->ref1 == it1->ref) && (fit->ref2 == it2->ref)) {
					fit->diff += it1->diff * it2->diff * factor;
					found = 1;
					break;
				}
				
			}

			if (found==0) {
				Magic2Deriv d;
				
				d.ref1 = it1->ref;
				d.ref2 = it2->ref;
			
				d.diff = it1->diff * it2->diff * factor;
				Stats[9]++;
				depends2.push_back(d);
			}

		}
	}
}



std::ostream &operator<<(std::ostream& os, const MagicDouble& a) {

	os << a.val << "[" << &a << "]" << a.depends << a.depends2;
	return os;

}


std::ostream &operator<<(std::ostream& os, const MagicDeriv& a) {

	os << "{" << a.ref << "->" << a.diff << "}";
	return os;

}

std::ostream &operator<<(std::ostream& os, const Magic2Deriv& a) {

	if (a.ref1==a.ref2) 
		os << "{{" << a.ref1 << "^2" << "->" << a.diff << "}}";
	else
		os << "{{" << a.ref1 << " " <<a.ref2 << "->" << a.diff << "}}";
	return os;

}

MagicDouble operator+(const MagicDouble& a, const MagicDouble& b) {

	MagicDouble r;
	r.val = a.val + b.val;
	r.depends_on(a.depends,1);
	r.depends_also_on(b.depends,1);

	r.depends2_on(a.depends2,1);
	r.depends2_also_on(b.depends2,1);
	
	return r;
}

MagicDouble& MagicDouble::operator+=(const MagicDouble& other) {
	
	val += other.val;
	depends_also_on(other.depends,1);
	depends2_also_on(other.depends2,1);
	
	return *this;
}

MagicDouble operator+(double a, const MagicDouble& b) {

	MagicDouble r;
	r.val = a + b.val;
	r.depends_on(b.depends,1);

	r.depends2_on(b.depends2,1);
	
	return r;
}

MagicDouble operator+(const MagicDouble& a, double b) {

	MagicDouble r;
	r.val = a.val + b;
	r.depends_on(a.depends,1);

	r.depends2_on(a.depends2,1);
	
	return r;
}

MagicDouble operator+(const MagicDouble& a) {

	MagicDouble r;
	r.val = a.val;
	r.depends_on(a.depends,1);

	r.depends2_on(a.depends2,1);
	
	return r;
}



MagicDouble operator-(const MagicDouble& a, const MagicDouble& b) {

	MagicDouble r;
	r.val = a.val - b.val;
	r.depends_on(a.depends,1);
	r.depends_also_on(b.depends,-1);

	r.depends2_on(a.depends2,1);
	r.depends2_also_on(b.depends2,-1);
	
	return r;
}

MagicDouble operator-(double a, const MagicDouble& b) {

	MagicDouble r;
	r.val = a - b.val;
	r.depends_on(b.depends,-1);

	r.depends2_on(b.depends2,-1);
	
	return r;
}

MagicDouble operator-(const MagicDouble& a, double b) {

	MagicDouble r;
	r.val = a.val - b;
	r.depends_on(a.depends,1);

	r.depends2_on(a.depends2,1);
	
	return r;
}

MagicDouble operator-(const MagicDouble& a) {

	MagicDouble r;
	r.val = -a.val;
	r.depends_on(a.depends,-1);

	r.depends2_on(a.depends2,-1);
	
	return r;
}

MagicDouble operator*(const MagicDouble& a, const MagicDouble& b) {

	MagicDouble r;
	r.val = a.val * b.val;
	r.depends_on(a.depends,b.val);
	r.depends_also_on(b.depends,a.val);

	r.depends2_on(a.depends2,b.val);           // x'' * y
	r.depends2_also_on(b.depends2,a.val);      // x * y''
	r.depends2_also_on(a.depends,b.depends,1); // x' * y'
	r.depends2_also_on(b.depends,a.depends,1); // y' * x'

	return r;
}

MagicDouble operator*(double a, const MagicDouble& b) {

	MagicDouble r;
	r.val = a * b.val;
	r.depends_on(b.depends,a);

	r.depends2_also_on(b.depends2,a);

	return r;
}

MagicDouble operator*(const MagicDouble& a, double b) {

	MagicDouble r;
	r.val = a.val * b;
	r.depends_on(a.depends,b);

	r.depends2_also_on(a.depends2,b);
	
	return r;
}

MagicDouble& MagicDouble::operator*=(double other) {

	val *= other;
	
	scale(other);
	scale2(other);
	
	return *this;
}

MagicDouble operator/(const MagicDouble& a, const MagicDouble& b) {

	MagicDouble r;
	r.val = a.val / b.val;
	r.depends_on(a.depends,1./b.val);
	r.depends_also_on(b.depends,-a.val/(b.val*b.val));

	r.depends2_on(a.depends2, 1/b.val);           // x'' 
	r.depends2_also_on(b.depends2, -a.val/(b.val*b.val) );      //  y''
	
	r.depends2_also_on(a.depends,b.depends, -1/(b.val*b.val)); // x' * y'
	r.depends2_also_on(b.depends,a.depends, -1/(b.val*b.val)); // x' * y'
	r.depends2_also_on(b.depends,b.depends, a.val*2/(b.val*b.val*b.val)); // y' * y'
	
	return r;
}

MagicDouble operator/(double a, const MagicDouble& b) {

	MagicDouble r;
	r.val = a / b.val;
	r.depends_on(b.depends,-a/(b.val*b.val));

	 
	r.depends2_on(b.depends2, -a/(b.val*b.val) );      //  y''
	r.depends2_also_on(b.depends,b.depends, a*2/(b.val*b.val*b.val)); // y' * y'

	
	return r;
}

MagicDouble operator/(const MagicDouble& a, double b) {

	MagicDouble r;
	r.val = a.val / b;
	r.depends_on(a.depends,1./b);

	r.depends2_on(a.depends2, 1/b);           // x'' 
	
	return r;
}

MagicDouble sin(const MagicDouble& a) {

	MagicDouble r;
	r.val = std::sin(a.val);
	r.depends_on(a.depends,cos(a.val));
	
	r.depends2_on(a.depends2,cos(a.val));
	r.depends2_also_on(a.depends,a.depends, -std::sin(a.val)); // y' * y'

	return r;
}

MagicDouble cos(const MagicDouble& a) {

	MagicDouble r;
	r.val = std::cos(a.val);
	r.depends_on(a.depends,-std::sin(a.val));
	
	r.depends2_on(a.depends2,-std::sin(a.val));
	r.depends2_also_on(a.depends,a.depends, -std::cos(a.val)); // y' * y'

	return r;
}

MagicDouble exp(const MagicDouble& a) {

	MagicDouble r;
	r.val = std::exp(a.val);
	r.depends_on(a.depends,std::exp(a.val));
	
	r.depends2_on(a.depends2,std::exp(a.val));
	r.depends2_also_on(a.depends,a.depends, std::exp(a.val)); // y' * y'

	return r;
}

}
