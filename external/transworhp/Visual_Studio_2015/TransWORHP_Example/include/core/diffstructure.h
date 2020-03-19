#pragma once

#include "../base/defines.h"
#include "../base/point.h"

#include <iostream>
#include <map>

#ifndef DllExport
#ifdef _WIN32
#define DllExport __declspec( dllexport )
#else
#define DllExport
#endif
#endif

#ifdef _MSC_VER
#pragma warning(disable : 4251)
#endif

namespace tw {

/** (dichtes) Abspeichern von Strukturinformationen und Ableitungswerten */
class DllExport DiffStructure {
	friend DllExport std::ostream &operator<<(std::ostream &os, const DiffStructure &d);

public:
	/**
	 * Constructor
	 */
	DiffStructure();

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

	/**
	 * explTW only! Set sparsity structure for non-multinodes
	 */
	void operator()(int i, int j, int dis, int n_param);


	int activeindex;

	bool ode_power(int n_ode, int n_ctrl);
	bool neben_power(DiffStructure &DS_ode, int n_ode);

	bool isValid() const;
	
	bool useStructure() const;
	void useStructure(bool use);


	// Fuer WORHP-Tool
	int GetEq() const;
	int GetDiff() const;
	int GetEntries() const;
	int CheckEntry(int i, int j) const;
	int check(int l, int j) const;
	double* find(int i, int j);
	double get(int i, int j) const;
	void finish();
	int odecheck(int j, int l, int block_index, int block_eq, double &val, double &mult, TWdiscretizationType type);
	
private:

	struct DiffEntry {
		int use;
		double value[3];
	};

	int n_eq, n_diff;
	bool adding;
	bool use_struct;
	std::map<Point<int>, DiffEntry> entrymap;
};

}
