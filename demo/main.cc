#include <chrono>
#include <iostream>
#include <malloc/_malloc.h>
#include <stdio.h>
#include <stdlib.h>

#include "odepack/dlsoda.h"

using namespace std;

struct PkData {
  double cl;
  double v;
  double ka;
};

int pk(int *neq, double *t, double *y, double *ydot, void *dat) {

  auto pk_data = *static_cast<PkData *>(dat);
  // auto pk_data = PkData{.cl = 0.134, .v = 7.81, .ka = 0.571};

  ydot[0] = -pk_data.ka * y[0];
  ydot[1] = pk_data.ka * y[0] - pk_data.cl / pk_data.v * y[1];
  return (0);
}

// int pk(int *neq, double *t, double *y, double *ydot,) {

//   // auto pk_data = *static_cast<PkData *>(dat);
//   auto pk_data = PkData{.cl = 0.134, .v = 7.81, .ka = 0.571};

//   ydot[0] = -pk_data.ka * y[0];
//   ydot[1] = pk_data.ka * y[0] - pk_data.cl / pk_data.v * y[1];
//   return (0);
// }

int test_lsoda() {

  int i, j, k, l, m, n, neq, itol, itask, istate, iopt, lrw, liw, jt, *iwork;
  double rtol, *atol, t, tout, *y, *yout, *rwork;

  neq = 2;
  iwork = (int *)malloc((neq + 20) * sizeof(int));
  y = (double *)malloc(neq * sizeof(double));
  yout = (double *)malloc(neq * sizeof(double));
  rwork = (double *)malloc((22 + neq * 16) * sizeof(double));
  atol = (double *)malloc(neq * sizeof(double));
  lrw = 22 + neq * 16;
  liw = 20 + neq;
  iopt = 0;
  jt = 2;
  t = 0.0E0;
  tout = 0.1;
  rtol = 1.0E-6;
  atol[0] = 1e-6;
  atol[1] = 1e-6;
  itol = 2;
  itask = 1;
  istate = 1;

  y[0] = 100.0;
  y[1] = 0.0;

  auto pk_data = PkData{.cl = 0.134, .v = 7.81, .ka = 0.571};

  for (i = 0; i < 12; i++) {
    dlsoda(pk, &neq, y, &t, &tout, &itol, &rtol, atol, &itask, &istate, &iopt,
           rwork, &lrw, iwork, &liw, 0, &jt, &pk_data);
    printf(" at t= %12.4e y= %14.6e %14.6e\n", t, y[0], y[1]);
    tout = tout + i;
  }

  free(y);
  free(iwork);
  free(yout);
  free(rwork);
  free(atol);

  return (0);
}

int main(void) {
  chrono::steady_clock::time_point begin, end;
  begin = chrono::steady_clock::now();
  double dt;

  int N = 1;

  // for (int i = 0; i < N; i++) {
  //   test_lsoda();
  // }
  // end = chrono::steady_clock::now();
  // dt = chrono::duration_cast<chrono::microseconds>(end - begin).count();
  // cout << "Time taken (us/per loop) = " << dt / N << endl;

  for (int i = 0; i < N; i++) {
    test_lsoda();
  }
  end = chrono::steady_clock::now();
  dt = chrono::duration_cast<chrono::microseconds>(end - begin).count();
  cout << "Time taken (us/per loop) = " << dt / N << endl;

  return (0);
}

/*
 The correct answer (up to certain precision):

 at t=   4.0000e-01 y=   9.851712e-01   3.386380e-05   1.479493e-02
 at t=   4.0000e+00 y=   9.055333e-01   2.240655e-05   9.444430e-02
 at t=   4.0000e+01 y=   7.158403e-01   9.186334e-06   2.841505e-01
 at t=   4.0000e+02 y=   4.505250e-01   3.222964e-06   5.494717e-01
 at t=   4.0000e+03 y=   1.831976e-01   8.941773e-07   8.168015e-01
 at t=   4.0000e+04 y=   3.898729e-02   1.621940e-07   9.610125e-01
 at t=   4.0000e+05 y=   4.936362e-03   1.984221e-08   9.950636e-01
 at t=   4.0000e+06 y=   5.161833e-04   2.065787e-09   9.994838e-01
 at t=   4.0000e+07 y=   5.179804e-05   2.072027e-10   9.999482e-01
 at t=   4.0000e+08 y=   5.283675e-06   2.113481e-11   9.999947e-01
 at t=   4.0000e+09 y=   4.658667e-07   1.863468e-12   9.999995e-01
 at t=   4.0000e+10 y=   1.431100e-08   5.724404e-14   1.000000e+00
 */
