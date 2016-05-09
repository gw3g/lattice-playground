#include <stdio.h>
#include "matrix.h"

/*
struct cuComplex {
  float Re; float Im;
  // constructors
  __device__ cuComplex()                             : Re(0.), Im(0.) {}
  __device__ cuComplex( float a, float b )           : Re( a), Im( b) {}
  // functions
  __device__ float     abs2( void )                  { return Re*Re + Im*Im; }
  __device__ cuComplex operator+(const cuComplex& a) { return cuComplex( Re+a.Re, Im+a.Im  ); }
  __device__ cuComplex operator-(const cuComplex& a) { return cuComplex( Re-a.Re, Im-a.Im  ); }
  __device__ cuComplex operator*(const cuComplex& a) { return cuComplex( Re*a.Re - Im*a.Im, 
                                                                         Im*a.Re + Re*a.Im ); }
  __device__ cuComplex operator*(const float a)      { return cuComplex( Re*a,  Im*a       ); }
};
*/

__device__ cuComplex conj(cuComplex a) {
  cuComplex ca(a.Re,-a.Im); return ca;
}
__device__ cuComplex inv(cuComplex a) {
  float c2 = a.abs2();
  cuComplex c(a.Re/c2,-a.Im/c2)    ; return c;
}


__global__ void view_m(cuComplex X[Nc][Nc]) {
  cuComplex tmp;
  for (int i=0; i<Nc; i++) {
    for (int j=0; j<Nc; j++) {
      tmp = X[i][j]; 
      printf("   (%.3f, %.3f) ", tmp.Re, tmp.Im );
    } printf("\n");
  }   printf("\n");
}

__device__ void copy_m(cuComplex A[Nc][Nc],  cuComplex B[Nc][Nc]) {
  for (int i=0; i<Nc; i++) 
    for (int j=0; j<Nc; j++) 
      B[i][j] = A[i][j];
}

__device__ void mul_m(cuComplex A[Nc][Nc],  cuComplex B[Nc][Nc], 
                                    cuComplex C[Nc][Nc]) {
  cuComplex Ct[Nc][Nc];
  for (int i=0; i<Nc; i++) 
    for (int j=0; j<Nc; j++) 
      for (int k=0; k<Nc; k++) Ct[i][j] = Ct[i][j] + A[i][k]*B[k][j];
  for (int i=0; i<Nc; i++) 
    for (int j=0; j<Nc; j++) 
      C[i][j] = Ct[i][j];
}

__device__ void add_m(cuComplex A[Nc][Nc],  cuComplex B[Nc][Nc], 
                                    cuComplex C[Nc][Nc]) {
  for (int i=0; i<Nc; i++) 
    for (int j=0; j<Nc; j++) C[i][j] = A[i][j]+B[i][j];
}

__device__ void sub_m(cuComplex A[Nc][Nc],  cuComplex B[Nc][Nc], 
                                    cuComplex C[Nc][Nc]) {
  for (int i=0; i<Nc; i++) 
    for (int j=0; j<Nc; j++) C[i][j] = A[i][j]-B[i][j];
}

__device__ void equ_m( cuComplex a,
                                    cuComplex A[Nc][Nc]) {
  cuComplex z(0.,0.);
  for (int i=0; i<Nc; i++) 
    for (int j=0; j<Nc; j++) A[i][j] = ((i==j) ? a : z);
}


__device__ void dag_m(cuComplex A[Nc][Nc]) {
  cuComplex temp;
  for (int i=0; i<Nc; i++)
    for (int j=0; j<=i; j++) {
      temp = A[i][j];
      A[i][j] = conj(A[j][i]);
      A[j][i] = conj(temp);
    }
}

__device__ cuComplex det_m(cuComplex A[Nc][Nc]) {
  cuComplex d;
  cuComplex zer[Nc][Nc];
  cuComplex Ap[Nc][Nc];
  add_m( zer, A, Ap);
  for (int i=0; i<Nc-1; i++)
    for (int j=i+1; j<Nc; j++)
      // subtract ~ row i from row j so that A[j][i] vanishes ...
      for (int k=i+1;k<Nc;k++) { Ap[j][k] = Ap[j][k] - ( Ap[j][i] * Ap[i][k]
          )*inv(Ap[i][i]); }
  d = Ap[0][0];
  for (int i=1; i<Nc; i++) d= d * Ap[i][i];
  return d;
}

__device__ cuComplex trace( cuComplex A[Nc][Nc] ) {
  cuComplex tr = A[0][0];
  for (int i=1; i<Nc; i++) tr = tr + A[i][i];
  return tr;
}

__device__ void gramschmidt(cuComplex A[Nc][Nc]) {
  cuComplex dot;
  float norm;
  for (int i=0; i<Nc; i++) {
    // normalise ith row
    /*norm = cabs( A[i][0] ); norm = norm * norm;*/
    norm = A[i][0].abs2();
    /*for (int j=1; j<Nc; j++) norm      = norm + cabs( A[i][j] )*cabs( A[i][j] );*/
    for (int j=1; j<Nc; j++) norm       += A[i][j].abs2();
    norm = 1./sqrt(norm);
    for (int j=0; j<Nc; j++) A[i][j]   = A[i][j]*norm;
    // orthogonalise the rest
    for (int k=i+1; k<Nc; k++) {
      dot = conj(A[i][0])*(A[k][0]);
      for (int j=1; j<Nc; j++) dot     = dot + conj(A[i][j])*(A[k][j]);
      for (int j=0; j<Nc; j++) A[k][j] = A[k][j] - dot*A[i][j];
    }
  }
}

__global__ void suN_m(cuComplex A[Nc][Nc]) {
  gramschmidt( A );
  cuComplex d; 
  /*int j, k;*/
  switch (Nc) { 
    // can do Nc=2,3 by hand
    case 2: 
      A[1][0] = conj(A[0][1])*(-1.);
      A[1][1] = conj(A[0][0])*(+1.);
      break;
    case 3:
      A[2][0] = conj( A[0][1]*A[1][2] - A[1][1]*A[0][2] );
      A[2][1] = conj( A[0][2]*A[1][0] - A[1][2]*A[0][0] );
      A[2][2] = conj( A[0][0]*A[1][1] - A[1][0]*A[0][1] );
      break; 
    default:
      d = det_m(A);
      for (int i=0; i<Nc; i++) A[0][i] = A[0][i] * conj(d);
  }
}

