#ifndef diffstructure_h
#define diffstructure_h

#include "../base/defines.h"
#include <iostream>
#include <map>

#include "../base/point.h"

class TransWorhp;


#ifdef WIN32
#ifndef MINGW
#define DllExport __declspec( dllexport )
#pragma warning (disable : 4251)
#else 
#define DllExport
#endif
#else 
#define DllExport
#endif


/** (dichtes) Abspeichern von Strukturinformationen und Ableitungswerten */
class DllExport DiffStructure {
	friend DllExport std::ostream &operator<<(std::ostream &os, const DiffStructure &d);
	friend class TransWorhp;
	friend class MagicTransWorhp;
	friend class ExplTransWorhp;
	friend class TWfolder;
	
public:
	/**
	 * Constructor
	 */
	DiffStructure();

	DiffStructure(const DiffStructure &other);
	/**
	 *
	 */
	~DiffStructure();

	/**
	 * Initialisierung
	 * @param eq Number of equations
	 * @param diff Number of dependency parameters
	 */
	void Init(int eq, int diff);
	
	/**
	 * Set sparsity structure (and set values)
	 * @param i equation index
	 * @param j dependency parameter index
	 */
	double& operator()(int i, int j);
	
	int activeindex;
	
	bool odepower(int n_ode, int n_ctrl);

	bool isValid();


	// Für WORHP-Tool
	int GetEq() {return n_eq;}
	int GetDiff() {return n_diff;}
	int GetEntries() {return (int) entrymap.size();}
	int CheckEntry(int i, int j) const {return check(i,j);}

private:
	
	double get(int i, int j);

	int adding;
	void finish();

	int odecheck(int j, int l, int block_index, int block_eq, double &val, double &mult, DiscretizationType_e type);
	int check(int l, int j) const;
	
	int n_eq, n_diff;

	struct DiffEntry {
		int use;
		double value[3];
	};

	double* find(int i, int j);

	std::map< Point<int>,  DiffEntry > entrymap;

	bool use_struct;
//	int *sum;
	
};

#endif
