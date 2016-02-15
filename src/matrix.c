/*
 * Author: greg jackson
 * Date: Feb 08 2016
 * Functions for complex Nc*Nc matrices. Intended for SU(Nc) calculations.
 */

#include "../include/matrix.h"

  void view_m(double complex  X[Nc][Nc]) {
  double complex tmp;
  for (int i=0; i<Nc; i++) {
    for (int j=0; j<Nc; j++) {
      tmp = X[i][j]; 
      printf("   (%.3f, %.3f) ", creal(tmp), cimag(tmp) );
    } printf("\n");
  }   printf("\n");
}

void mul_m(double complex A[Nc][Nc],  double complex B[Nc][Nc], 
                                    double complex C[Nc][Nc]) {
  for (int i=0; i<Nc; i++) 
    for (int j=0; j<Nc; j++) 
      for (int k=0; k<Nc; k++) C[i][j] += A[i][k]*B[k][j];
}

void add_m(double complex A[Nc][Nc],  double complex B[Nc][Nc], 
                                    double complex C[Nc][Nc]) {
  for (int i=0; i<Nc; i++) 
    for (int j=0; j<Nc; j++) C[i][j] = A[i][j]+B[i][j];
}

void sub_m(double complex A[Nc][Nc],  double complex B[Nc][Nc], 
                                    double complex C[Nc][Nc]) {
  for (int i=0; i<Nc; i++) 
    for (int j=0; j<Nc; j++) C[i][j] = A[i][j]-B[i][j];
}

void equ_m( double complex a,
                                    double complex A[Nc][Nc]) {
  for (int i=0; i<Nc; i++) 
    for (int j=0; j<Nc; j++) A[i][j] = ((i==j) ? a : .0);
}

void conj_m(double complex A[Nc][Nc]) {
  for (int i=0; i<Nc; i++)
    for (int j=0; j<Nc; j++) A[i][j] = conj(A[i][j]);
}

double complex det_m(double complex A[Nc][Nc]) {
  double complex d;
  double complex zer[Nc][Nc] = { {.0} }, Ap[Nc][Nc];
  add_m( zer, A, Ap);
  for (int i=0; i<Nc-1; i++)
    for (int j=i+1; j<Nc; j++)
      // subtract ~ row i from row j so that A[j][i] vanishes ...
      for (int k=i+1;k<Nc;k++) { Ap[j][k] -= Ap[j][i]*Ap[i][k]/Ap[i][i]; }
  d = Ap[0][0];
  for (int i=1; i<Nc; i++) d*=Ap[i][i];
  return d;
}

void gramschmidt(double complex A[Nc][Nc]) {
  double complex dot, norm;
  for (int i=0; i<Nc; i++) {
    // normalise ith row
    norm = cabs( A[i][0] ); norm*=norm;
    for (int j=1; j<Nc; j++) norm      += cabs( A[i][j] )*cabs( A[i][j] );
    norm = 1./csqrt(norm);
    for (int j=0; j<Nc; j++) A[i][j]   *= norm;
    // orthogonalise the rest
    for (int k=i+1; k<Nc; k++) {
      dot = conj(A[i][0])*(A[k][0]);
      for (int j=1; j<Nc; j++) dot     += conj(A[i][j])*(A[k][j]);
      for (int j=0; j<Nc; j++) A[k][j] -= dot*A[i][j];
    }
  }
}

void suN_m(double complex A[Nc][Nc]) {
  gramschmidt( A );
  double complex d; 
  /*int j, k;*/
  switch (Nc) { 
    // can do Nc=2,3 by hand
    /*
    case 2: 
      A[1][0] = -conj(A[0][1]);
      A[1][1] = +conj(A[0][0]);
      break;
    case 3:
      A[2][0] = conj( A[0][1]*A[1][2] - A[1][1]*A[0][2] );
      A[2][1] = conj( A[0][2]*A[1][0] - A[1][2]*A[0][0] );
      A[2][2] = conj( A[0][0]*A[1][1] - A[1][0]*A[0][1] );
      break; 
      */
    default:
      d = det_m(A);
      for (int i=0; i<Nc; i++) A[0][i] *= conj(d);
  }
}

