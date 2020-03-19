#ifndef CPP_CS_H
#define CPP_CS_H

extern "C" {
#include "macros.h"
#include "C_cs.h"

  DLL_PUBLIC void SortWorhpMatrix(WorhpMatrix *WM);
  void Reorder( mat_int *v, mat_int *order, size_t size );
}

/*
 * Define struct for sort, this way
 * we may provide an object of WorhpMatrixSort
 * to sort and therefore obtain values inside of WM
 *
 * @AUTHOR: sge
 */
struct WorhpMatrixSort
{
  WorhpMatrix* WM;
  WorhpMatrixSort(WorhpMatrix* p) : WM(p) {};

  bool operator() ( const int &i, const int &j )
  {
    if ( WM->kind != Matrix_Kind_LowTri ) {
      if (WM->col[i-1] == WM->col[j-1]) {
	return WM->row[i-1] <= WM->row[j-1];
      } else {
	return WM->col[i-1] <= WM->col[j-1];
      }
    } else {
      if (WM->col[i-1] != WM->row[i-1] && WM->col[j-1] != WM->row[j-1]) {
	if (WM->col[i-1] == WM->col[j-1]) return WM->row[i-1] <= WM->row[j-1];
	else return WM->col[i-1] <= WM->col[j-1];	
      }
      else if (WM->col[i-1] == WM->row[i-1] && WM->col[j-1] == WM->row[j-1]) {
	return WM->col[i-1] < WM->col[j-1];
      }
      else if (WM->col[i-1] == WM->row[i-1] && WM->col[j-1] != WM->row[j-1]) {
	return 0;
      }
      else /*if (col!=row && other.col==other.row)*/ {
	return 1;
      }    
    }
  }
};
#endif
