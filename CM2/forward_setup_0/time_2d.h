#ifndef _TIME_2D_
#define _TIME_2D_


  #include <stdlib.h>
  int time_2d(double *HS, double *T, int NX, int NY, double XS, double YS, double EPS_INIT, int MESSAGES);
  void scaling(double* hsbuf, double* realxs, double* realys, double* S, double h, int size_S);
  void connect();
#endif /* _TIME_2D_ */
