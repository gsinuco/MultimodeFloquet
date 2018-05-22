#include "constants.h"

extern"C" {
 // void  rwaresonantfrequencies_(double *DC_FIELD,double *RF_FIELD,double *MW_FIELD,double *RESONAN_FREQUENCIES,double *MW_COUPLINGS,int *INFO);
 void        rwaresonantfrequencies_(int *t_,double *DC_FIELD,double *RF_FIELD,double *MW_FIELD,double *RESONAN_FREQUENCIES,double *MW_COUPLINGS,int *INFO);
 void    floquetresonantfrequencies_(int *t_,double *DC_FIELD,double *RF_FIELD,double *MW_FIELD,double *RESONAN_FREQUENCIES,double *MW_COUPLINGS,int *INFO);
 void  multimodefloquetimeevolution_(int *t_,double *DC_FIELD,double *RF_FIELD,double *MW_FIELD,double *T1,double *T2,double *Pup,int *INFO);
}


extern"C" {
  void floquetinit();
  void lapack_fulleigenvalues();
  void dressedbasis();
  void sethamiltoniancomponents();
  void multimodefloquetmatrix();
  void multimodefloquetmatrix_sp();
  void mklsparse_fulleigenvalues();
  void multimodetimeevolutionoperator();
  void multimodefloquettransformation();
}
