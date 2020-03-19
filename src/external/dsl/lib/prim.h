#ifndef PRIM_H
#define PRIM_H

#include <vector>

/**
 * This header defines the prim struct and contains all information about the prim these are the states,
 * position and heading, as well as the dynamic properties and the controls of the trim. Some additional
 * information about the primitive.
 */ 

struct prim {

	bool is_trim;
	int id;
	int n_t_steps;
	double costs;
	int from;
	int to;
	std::vector<double> t;
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> psi;
	std::vector<double> u;
	std::vector<double> v;
	std::vector<double> r;
	//neu
	std::vector<double> delta;
	std::vector<double> n;
	std::vector<double> energy;
	std::vector<double> d_delta;
	std::vector<double> d_n;
	//
	std::vector<double> d;

};

#endif