#include "TWstrings.h"

std::string TWstrings::string_hessianstructure[10] = {
	"(structure=0) Using Diagonal Matrix",
	"(structure=1) Using Full Matrix",
	"(structure=2) Using oversized ode blocks",
	"(structure=3) Using ode blocks ",
	"(structure=4) Using analytic derivative TOO SLOW!!!",
	"(structure=5) Using analytic derivative+++ ",
	""
};

std::string TWstrings::string_hessianvalues[10] = {
	"(values=0) Calculating 2nd derivative of F + mu*G",
	"(values=1) Calculating 1st derivative of DF + mu*DG",
	"(values=2) Calculating 1st derivative of DF + mu*DG ... faster",
	"(values=3) Calculating analytic derivative of DF + mu*DG",
	"(values=4) Calculating analytic derivative of DF + mu*DG 2nd version",
	"",
};

// abhaengig von DiscretizationType
std::string TWstrings::string_discretization[20] = {
	"Using Trapez discretization",
	"Using HermiteSimpson discretization",
	"Using Explicit Euler discretization",
	"Using Lobatto discretization",
	"",
	"",
	"",
	"",
	"",
	"",
	"Using Single/Multiple Shooting", // 10
	"Using Pseudospectral Method",
	"",
	"",
};
