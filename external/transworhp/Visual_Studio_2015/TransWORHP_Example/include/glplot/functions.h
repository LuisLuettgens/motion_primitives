#pragma once

/**
 *  Use this definition to implement your own data drawing routine, see also
 *  createuserwindow_().
 *  @param handle Window handle.
 *  @param frame Current frame to be drawn.
 *  @param len Number of available data sets.
 *  @param ndgl Dimension of x data set.
 *  @param nsteuer Dimension of u data set.
 *  @param t Time vector (of size len).
 *  @param x ODE vector (of size len*ndgl).
 *  @param u Control vector (of size len*nsteuer).
 *  @param bord Set this double[4]-Vector to provide min/max values for axes.
 */
using Funktionenzeiger = void (*) (int &handle, int &frame,
                                  int &len, int &ndgl, int &nsteuer, double *t, double *x,
                                  double *u, double *bord);
using FunktionenzeigerI = void (*) (int &handle, int &frame,
                                   int &len, int &ndgl, int &nsteuer, double *t, double *x,
                                   double *u, double *bord, int &index);
using FunktionenzeigerU = void (*) (int &handle, int &frame,
                                   int &len, int &ndgl, int &nsteuer, int &nunbe, double *t, double *x,
                                   double *u, double *unknown, double *bord);
/**
 *  Use this definition to draw more complex data sets by combining values from
 *  t, x and u, see also createdatawindow_().
 *  @param len Number of available data sets.
 *  @param ndgl Dimension of x data set.
 *  @param nsteuer Dimension of u data set.
 *  @param t Time vector (of size len).
 *  @param x ODE vector (of size len*ndgl).
 *  @param u Control vector (of size len*nsteuer).
 *  @param i Data at position i (out of len) should be returned.
 *  @param index Specify type of value.
 *  @return Calculated data of type index at position i.
 */
using Funktionenzeiger2 = double (*) (int &len, int &ndgl, int &nsteuer,
                                     double *t,double *x, double *u, int &i, int &index);
/** @} */
using FunktionenzeigerC = void (*) (int &button, double &xpos, double &ypos,
                                   int &len, int &ndgl, int &nsteuer, double *t, double *x, double *u);

using timefuncf = double (*) (int &, double *, double *, double *, int &, int &, int &);
