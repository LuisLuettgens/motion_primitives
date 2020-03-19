#ifndef course_correction_H
#define course_correction_H

#include "shipviewer_base_Nships.h" 
#include "geometric_constr.h"

using namespace std;

void korrektur_startposition(vector<double> pertubation);

void korrektur_endposition(vector<double> shift);

void korrektur(vector<double> pertubation, vector<double> shift);

#endif