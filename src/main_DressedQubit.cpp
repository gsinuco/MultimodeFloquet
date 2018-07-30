#include <iostream>
#include <complex>
#include <stdio.h>
#include <math.h>
#include <string.h>

using namespace std;
typedef std::complex<double> dcmplx;

#include "MultimodeFloquet.h"

extern "C" int h_floquet_size;

int main(){

  atom_c id;
  int info;
  int jtotal;
  //  char *name = new char [6];  //[5] = "qubit";
  char name[]     = "qubit";
  char manifold[] = "U";
  char op[]="N";


  int r,m,l,i,j;
  int d_bare,total_frequencies,sp;

  double t1,t2,e_l,e_r;

  info   = 0;
  jtotal = 2;
  floquetinit_c(name,manifold,&jtotal,&id,&info);

  d_bare = id.d_bare;

  dcmplx * U_AUX = new dcmplx [d_bare*d_bare];
  double * P = new double [d_bare*d_bare];

  int nm = 3;
  int * modes_num = new int [nm];

  modes_num[0] = 1;
  modes_num[1] = 1;
  modes_num[2] = 1;
  
  total_frequencies = 0;
  for(r=0;r<nm;r++){
    total_frequencies += modes_num[r];
  }
  
  mode_c * fields = new mode_c [total_frequencies];
  
  
  fields[0].x     = 0.0;
  fields[0].y     = 0.0;
  fields[0].z     = 1.0;
  fields[0].phi_x = 0.0;
  fields[0].phi_y = 0.0;
  fields[0].phi_z = 0.0;
  fields[0].omega = 0.0;
  fields[0].N_Floquet = 0;

  fields[1].x     = 0.125;
  fields[1].y     = 0.0;
  fields[1].z     = 0.0;
  fields[1].phi_x = 0.0;
  fields[1].phi_y = 0.0;
  fields[1].phi_z = 0.0;
  fields[1].omega = 1.0;
  fields[1].N_Floquet = 5;

  fields[2].x     = 0.125*fields[1].x/2.0;
  fields[2].y     = 0.125*fields[1].x/2.0;
  fields[2].z     = 0.0;
  fields[2].phi_x = 0.0;
  fields[2].phi_y = 0.0;
  fields[2].phi_z = 0.0;
  fields[2].omega = real(fields[1].x)/2.0;
  fields[2].N_Floquet = 3;

  printf("%i %i \n",d_bare,total_frequencies);

  sethamiltoniancomponents_c_(&id,&nm,&total_frequencies,modes_num,fields,&info);
         
  // =================================================================================
  // ==== DEFINITION OF THE DRESSING FIELDS AND DRESSED BASIS
  // =================================================================================

  int dressingfields = 2;
  int * dressingfields_indices =  new int [dressingfields];
  
  dressingfields_indices[0] = 0;
  dressingfields_indices[1] = 1;
  
  int dressingfloquetdimension = d_bare;  
  for(m=0;m<dressingfields;m++){
      dressingfloquetdimension = dressingfloquetdimension*(2*fields[dressingfields_indices[m]].N_Floquet + 1);
  }
  dcmplx * U_FD = new dcmplx [dressingfloquetdimension*dressingfloquetdimension];
  double * e_dressed = new double [dressingfloquetdimension];
  
  dressedbasis_subset_c_(&id,&dressingfloquetdimension,&dressingfields,&nm,dressingfields_indices,modes_num,fields, U_FD, e_dressed,&info);

  
  int index0 = d_bare*fields[1].N_Floquet;
  
  int nm_ = dressingfields;
  int modes_num_[nm_];
  
  int total_frequencies_ = 0;
  for(r=0;r<nm_;r++){
      modes_num_[r] = modes_num[dressingfields_indices[r]];
      total_frequencies_ += modes_num_[r];
  }
  
  mode_c * fields_ = new mode_c [total_frequencies_];
  
  
  int field_index  = 0;
  int field_offset = 0;
  
  for(r=0;r<dressingfields;r++){
      field_offset = 0;
      for(l=0;l<dressingfields_indices[r];l++){
        field_offset += modes_num[l];
      }
      for(m=0;m<modes_num_[r];m++){
          fields_[r] = fields[field_offset + m];          
      }
      //printf("%i\n",field_offset);
  }
  
  
  //==== ALLOCATE ARRAYS FOR THE MICROMOTION OPERATORS IN THE EXTENDED AND ORIGINAL HILBERT SPACES.
  
  dcmplx * U_F1     = new dcmplx [d_bare*dressingfloquetdimension];
  dcmplx * U_F2     = new dcmplx [d_bare*dressingfloquetdimension];
  dcmplx * U_F1_red = new dcmplx [d_bare*d_bare];
  dcmplx * U_F2_red = new dcmplx [d_bare*d_bare];
  
  // ! ========= FIND THE MULTIMODE FLOQUET SPECTRUM 

  for(r=0;r<1;r++){
    fields[2].omega = real(fields[1].x)/4.0;// + r*fields[1].x/128;     
    sethamiltoniancomponents_c_(&id,&nm,&total_frequencies,modes_num,fields,&info); // every time a field parameter is modified, we should run this function
    //!--- FIND THE MULTIMODE FLOQUET SPECTRUM 
    
    multimodefloquetmatrix_c_(&id,&nm,&total_frequencies,modes_num,fields,&info); // in this function we calculate the dimension of the multimode floquet hilbert space, 
                                                                                  // which is the value of the global variable h_floquet_size    
    double * e_floquet = new double [h_floquet_size];
    dcmplx * U_F =  new dcmplx [h_floquet_size*h_floquet_size];    
    lapack_fulleigenvalues_c_(U_F,&h_floquet_size,e_floquet,&info);// here, the diagonalization is done with the internal (Fortran) Hamiltonian (H_FLOQUET)
                                                                   // U_F is the transformation that diagonalise the Hamiltonian
                                                                   // On the Fortran side, H_FLOQUET is deallocated after diagonalization. This is needed since
                                                                   // because the floquet hamiltonian is allocated again when running the subroutine
                                                                   // MULTIMODEFLOQUETMATRIX 
    //for(r=0;r<h_floquet_size;r++) printf("%15.5f  ",e_floquet[r]);
    //printf("\n");

    // ======= EVALUATE TIME-EVOLUTION OPERATOR IN THE BARE BASIS       
    t1 = 0.0;
    t2 = 1.0;
    for(m=0;m<1;m++){
        t2 = m*16.0*100/32.0;
        multimodetimeevolutionoperator_c_(&h_floquet_size,&nm,modes_num,U_F,e_floquet,&d_bare,fields,&t1,&t2,U_F,&info);
    //!=================================================================================
    //!== TRANSFORM THE TIME-EVOLUTION OPERATOR TO THE DRESSED BASIS
    //!=================================================================================
    //       
    //!== BUILD THE TIME-DEPENDENT TRANSFORMATION BETWEEN THE BARE AND THE RF DRESSED BASIS: U_F1
    //       
        info =0  ;
        multimodefloquettransformation_c_(&dressingfloquetdimension,&nm_,modes_num_,U_FD,e_dressed,&d_bare,fields_,&t1,U_F1,&info); 
        multimodefloquettransformation_c_(&dressingfloquetdimension,&nm_,modes_num_,U_FD,e_dressed,&d_bare,fields_,&t2,U_F2,&info); 
        /*for(i=0;i<d_bare*dressingfloquetdimension;i++){
            printf("%i %f \n ",i,abs(U_F1[i]));
        }*/
        
    //! ====== SINGLE OUT ONE BARE SUBSPACE
        
        index0 = d_bare*d_bare*(((dressingfloquetdimension/d_bare) - 1)/2);
        for(i=0;i<d_bare;i++){
            for(j=0;j<d_bare;j++){
            U_F1_red[i*d_bare+j] = U_F1[index0 + i*d_bare + j];
            U_F2_red[i*d_bare+j] = U_F2[index0 + i*d_bare + j];
            printf("%i %f \n",i*d_bare+j,abs(U_F2_red[i*d_bare+j]));
            }
        }
        printf("\n\n");
        multimodemicromotion_c_(&id,&dressingfloquetdimension,&nm_,modes_num_,U_FD,e_dressed,&d_bare,fields_,&t1,U_F1_red,&info); 
        for(i=0;i<d_bare;i++){
            for(j=0;j<d_bare;j++){
            printf("%i %f \n",i*d_bare+j,abs(U_F2_red[i*d_bare+j]));
            }
        }
 
    //! ====== CALCULATE THE TIME-EVOLUTION OPERATOR IN THE DRESSED BASIS USING THE PREVIOUSLY CALCULATED IN THE BARE BASIS
        matmul_c("N",U_F,&d_bare,&d_bare,U_F1_red,&d_bare,&d_bare,U_AUX,&info);
        matmul_c("TC",U_F2_red,&d_bare,&d_bare,U_AUX,&d_bare,&d_bare,U_F,&info);
        for(i=0;i<d_bare*d_bare;i++){
            P[i] = abs(U_F[i])*abs(U_F[i]);
        }
        
        write_matrix_c_(P,&d_bare);        
    }
    delete(e_floquet);
    delete(U_F);
  }

  return 0;
}


