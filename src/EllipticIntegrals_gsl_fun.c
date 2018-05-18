#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_sf_legendre.h>
/*
extern struct{
  double k;
  int l;
}EllipticArg_;
*/
void completek_(double *k, double *K){
  gsl_mode_t mode;
  mode = GSL_PREC_DOUBLE; 
  *K = gsl_sf_ellint_Kcomp(*k, mode);
}

void completee_(double *k, double *K){
  gsl_mode_t mode;
  mode = GSL_PREC_DOUBLE; 
  *K = gsl_sf_ellint_Ecomp(*k, mode);
}
/*
void legendre_pol_(double *k, int l, double *Pl){

  gsl_mode_t mode;
  mode = GSL_PREC_DOUBLE;
  *Pl =gsl_sf_legendre_Pl (l, *k);
  }*/
