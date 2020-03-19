#include "diffstructure.h"
#include <iomanip>
#include "worhp_info.h"
#include "conversion.h"

using namespace std;

namespace tw {

/** fixer Rueckgabewert */
double tmpdiff = 0;

DiffStructure::DiffStructure() :
	activeindex(0),
	n_eq(0), n_diff(0),
	adding(true),
	use_struct(true) {
	
}

void DiffStructure::Init(int eq, int diff) {
	
	n_eq = eq;
	n_diff = diff;
	
	entrymap.clear();
	adding = true;
}

int DiffStructure::odecheck(int j, int l, int block_index, int block_eq, double &val, double &mult, TWdiscretizationType type) {

	if (j>=n_eq)
		throw string("Bad Dim1 in Check:") + std::to_string(j) + ">=" + std::to_string(n_eq);
	if (l>=n_diff)
		throw string("Bad Dim2 in Check:") + std::to_string(l) + ">=" + std::to_string(n_diff);

	int pattern = 0;
	int eq = 0;
	val = 0;
	mult = 0;

	int data_ = 0;
	auto it = entrymap.find(Point<int>(j,l));
	if (it!=entrymap.end()) data_ = it->second.use;


	if (type == TWdiscretizationType::Trapez) {

		if (block_index==0) {
			if (j==l) {
				val=-1;
				eq = 1;
			}
			pattern = data_;
			mult=-.5;
		}
		if (block_index==1) {
			if (j==l) {
				val=1;
				eq = 1;
			}
			pattern = data_;
			mult=-.5;
		}


		//if (pattern) val=M_PI;
	}
	else if (type == TWdiscretizationType::HermiteSimpson) {

		if (block_eq==0) {
			if (block_index==0) {
				if (j==l) {
					val=-.5;
					eq = 1;
				}
				pattern = data_;
				mult=-.125;
			}
			if (block_index==1) {
				if (j==l) {
					val=1;
					eq = 1;
				}
				mult=0.;
				// pattern = data[j][l];
			}
			if (block_index==2) {
				if (j==l) {
					val=-.5;
					eq = 1;
				}
				pattern = data_;
				mult=.125;
			}
		}
		if (block_eq==1) {
			if (block_index==0) {
				if (j==l) {
					val=-1;
					eq = 1;
				}
				pattern = data_;
				mult=-1./6.;
			}
			if (block_index==1) {
				if (j==l) {}
				pattern = data_;
				mult=-4./6.;
			}
			if (block_index==2) {
				if (j==l) {
					val=1;
					eq = 1;
				}
				pattern = data_;
				mult=-1./6.;
			}
		}
		//if (pattern) val=M_PI;

	}
	else if (type == TWdiscretizationType::Lobatto) {
	
		if (block_eq==0) {
			if (block_index==0) {
				if (j==l) {
					val=-1./6+(5./6)*sqrt(5.);
					eq = 1;
				}
				pattern = data_;
				mult=-1./20.+(43./420.)*sqrt(5.);
			}
			if (block_index==1) {
				if (j==l) {
					val=-31./21;
					eq = 1;
				}
				mult=0.;
				// pattern = data[j][l];
			}
			if (block_index==2) {
				if (j==l) {
					val=23./14-(5./6)*sqrt(5.);
					eq = 1;
				}
				mult=0.;
				// pattern = data[j][l];
			}
			if (block_index==3) {
				if (j==l) {
					val=0;
					eq = 1;
				}
				pattern = data_;
				mult=(2/105.)*sqrt(5.)-1/30.;
			}
		}
		if (block_eq==1) {
			if (block_index==0) {
				if (j==l) {
					val=13/6.+(5/6.)*sqrt(5.);
					eq = 1;
				}
				pattern = data_;
				mult=7/30.+(3/35.)*sqrt(5.);
			}
			if (block_index==1) {
				if (j==l) {
					val=-23/14.-(5/6.)*sqrt(5.);
					eq = 1;
				}
				mult=0.;
				// pattern = data[j][l];
			}
			if (block_index==2) {
				if (j==l) {
					val=-11/21.;
					eq = 1;
				}
				mult=0.;
				// pattern = data[j][l];
			}
			if (block_index==3) {
				if (j==l) {
					val=0;
					eq = 1;
				}
				pattern = data_;
				mult=1/60.+(1/420.)*sqrt(5.);
			}
		}
		if (block_eq==1) {
			if (block_index==0) {
				if (j==l) {
					val=-1;
					eq = 1;
				}
				pattern = data_;
				mult=-1./12.;
			}
			if (block_index==1) {
				if (j==l) {}
				pattern = data_;
				mult=-5./12.;
			}
			if (block_index==2) {
				if (j==l) {}
				pattern = data_;
				mult=-5./12.;
			}
			if (block_index==3) {
				if (j==l) {
					val=1;
					eq = 1;
				}
				pattern = data_;
				mult=-1./12.;
			}
		}
		//if (pattern) val=M_PI;

	}
	else if (type == TWdiscretizationType::Euler) {

		if (block_index==0) {
			if (j==l) {
				val=-1;
				eq = 1;
			}
			pattern = data_;
			mult=-1;
		}
		if (block_index==1) {
			if (j==l) {
				val= 1;
				eq = 1;
			}
			//pattern = data_;
			mult=0;
		}

		//if (pattern) val=M_PI;
	}
	else if (type == TWdiscretizationType::MultipleShooting) {
		if (block_index==0) {
			if (j==l) {
				val= 1;
				eq = 1;
			}
			pattern = data_;
			mult=-1;
		}
		
	}


	if (use_struct && pattern==0) mult=0;


	return (eq || pattern || !use_struct);
}


int DiffStructure::check(int j, int l) const {

	if (j>=n_eq)
		throw string("Bad Dim1 in Check:") + std::to_string(j) + ">=" + std::to_string(n_eq);
	if (l>=n_diff)
		throw string("Bad Dim2 in Check:") + std::to_string(l) + ">=" + std::to_string(n_diff);

	int data_ = 0;
	auto it = entrymap.find(Point<int>(j,l));
	if (it!=entrymap.end()) data_ = it->second.use;

	return (!use_struct || data_);
}

bool DiffStructure::isValid() const {

	if (n_diff==0 || n_eq==0) {
		return false;
	} else {
		return true;
	}
}


bool DiffStructure::useStructure() const {
	return use_struct;
}


void DiffStructure::useStructure(bool use) {
	use_struct = use;
}


double DiffStructure::get(int i, int j) const {

	std::map< Point<int>,  DiffEntry >::const_iterator it = entrymap.find(Point<int>(i,j));
	if (it!=entrymap.end()) {
		return it->second.value[activeindex];
	} else {
		return 0;
	}
}


double& DiffStructure::operator()(int i, int j) {
	
	auto it = entrymap.find(Point<int>(i,j));
	if (it!=entrymap.end()) {
		return it->second.value[activeindex];
	} else {
		if (adding) {
			
			Point<int> k(i,j);
			DiffEntry de;
			de.use = 1;
			entrymap.insert ( pair<Point<int>, DiffEntry>(k,de) );
			it = entrymap.find(Point<int>(i,j));
			
			return it->second.value[activeindex];
		} else {
			cout << "Diffstructur: Index ("<<i<<","<<j<<") not found in structure " << this->n_eq << " " << this->n_diff << endl;
		}
		
		return tmpdiff;
	}
}

void DiffStructure::operator()(int i, int j, int dis, int n_param) {
	
	if (adding) {
		
		for (int k = j; k < j+dis; ++k) {
			Point<int> p(i,k);
			DiffEntry de;
			de.use = 1;
			entrymap.insert ( pair<Point<int>, DiffEntry>(p,de) );
		}
		
		for (int k = n_diff-n_param; k < n_diff; ++k) {
			Point<int> p(i,k);
			DiffEntry de;
			de.use = 1;
			entrymap.insert ( pair<Point<int>, DiffEntry>(p,de) );
		}
	}
}


double* DiffStructure::find(int i, int j) {

	auto it = entrymap.find(Point<int>(i,j));
	if (it!=entrymap.end()) {
		return it->second.value;
	} else {
		return nullptr;
	}
}


void DiffStructure::finish() {

	if (!use_struct) {
		for (int i = 0; i < n_diff; i++) {
			for (int j = 0; j < n_eq; j++) {
				(*this)(j,i);
			}
		}
	}

	adding = false;
}


bool DiffStructure::ode_power(int n_ode, int n_ctrl) {
	
	//cout << "Calc. power of ode structure" << endl;
	bool changed = false;
	
	for (int j=0;j<n_ode;j++) {
		auto itl = entrymap.find(Point<int>(j,j));
		if ( itl==entrymap.end() ) {
			Point<int> kk(j,j);
			DiffEntry de;
			de.use = 1;
			entrymap.insert ( pair<Point<int>, DiffEntry>(kk,de) );
			changed = true;
		}
	}
	
	for (int i=0;i<n_ode;i++) {
		for (int k=0;k<n_ode;k++) {
			
			auto it = entrymap.find(Point<int>(i,k));
			
			if (it==entrymap.end()) {
				
				for (int j=0;j<n_ode;j++) {
					auto itl = entrymap.find(Point<int>(i,j));
					auto itr = entrymap.find(Point<int>(j,k));
					if ( itl!=entrymap.end() && itr!=entrymap.end()) {
						
						Point<int> kk(i,k);
						DiffEntry de;
						de.use = 1;
						entrymap.insert ( pair<Point<int>, DiffEntry>(kk,de) );
						changed = true;
						break;
					}
				}
			}
		}
	}
	
	for (int i=0;i<n_ode;i++) {
		for (int k=0;k<n_ctrl;k++) {
			
			auto it = entrymap.find(Point<int>(i,n_ode + k));
			
			if (it==entrymap.end()) {
				
				for (int j=0;j<n_ode;j++) {
					auto itl = entrymap.find(Point<int>(i,j));
					auto itr = entrymap.find(Point<int>(j,n_ode + k));
					
					if (itl!=entrymap.end() && itr!=entrymap.end()) {
						
						Point<int> kk(i,n_ode + k);
						DiffEntry de;
						de.use = 1;
						entrymap.insert ( pair<Point<int>, DiffEntry>(kk,de) );
						changed = true;
						break;
					}
				}
			}
			
		}
	}
	
	return changed;
}

bool DiffStructure::neben_power(DiffStructure &DS_ode, int n_ode) {
	
	//cout << "Calc. structur of neben" << endl;
	bool changed = false;
	
	for (int i = 0; i < n_eq; ++i) {
		for (int k = 0; k < n_ode; ++k) {
			
			if (entrymap.find(Point<int>(i,k)) != entrymap.end()) {
				
				for (int x = 0; x < n_diff; ++x) {
					
					if (DS_ode.check(k,x)) {
						
						Point<int> kk(i,x);
						DiffEntry de;
						de.use = 1;
						entrymap.insert ( pair<Point<int>, DiffEntry>(kk,de) );
						changed = true;
					}
					
				}
			}
		}
	}
	
	return changed;
}


ostream &operator<<(ostream &os, const DiffStructure &d) {

	if (d.use_struct) {
		for (int i=0;i<d.n_eq;i++) {
			os << "   ";
			for (int j=0;j<d.n_diff;j++) {

				int data_ = -1;

				auto it = d.entrymap.find(Point<int>(i,j));
				if (it!=d.entrymap.end()) data_ = it->second.use;

				if (data_==-1) os << setw(2) << ".";
				else if (data_==0) os << setw(2) << "0";
				else os << setw(2) << "#";
			}
			os << endl;
		}
	}
	return os;
}

int DiffStructure::GetEq() const {
	return n_eq;
}

int DiffStructure::GetDiff() const {
	return n_diff;
}

int DiffStructure::GetEntries() const {
	return static_cast<int>(entrymap.size());
}

int DiffStructure::CheckEntry(int i, int j) const {
	return check(i,j);
}

}
