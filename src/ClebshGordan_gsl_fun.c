#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_coupling.h>

extern struct{
 double  l1,l2,l3,m1,m2,m3;
}jSymbols_;


void cb_fun_(int *l1,int *l2,int *l3,  int *m1, int *m2, int *m3, double *beta){  
  *beta = sqrt(2.0*0.5*(*l3) + 1.0)*pow(-1.0,0.5*(*l1 - *l2 + *m3))*gsl_sf_coupling_3j(((*l1)),((*l2)),((*l3)),((*m1)),((*m2)),(-(*m3)));  
  //  *beta = gsl_sf_coupling_3j(((*l1)),((*l2)),((*l3)),((*m1)),((*m2)),(-(*m3)));  
}

