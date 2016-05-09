#include <stdio.h>
#include <math.h>
#pragma once
#include <complex.h>

/*-----------------------------------------------------------------------------------------------*/

#define Nc 2

struct cuComplex {
  float Re; float Im;
  // constructors
  __device__ cuComplex()                             : Re(0.), Im(0.) {}
  __device__ cuComplex( float a, float b )           : Re( a), Im( b) {}
  // functions
  __device__ float     abs2( void )                  { return Re*Re + Im*Im; }
  __device__ cuComplex operator=(const float a)      { return cuComplex( a, 0. ); }
  __device__ cuComplex operator+(const cuComplex& a) { return cuComplex( Re+a.Re, Im+a.Im  ); }
  __device__ cuComplex operator-(const cuComplex& a) { return cuComplex( Re-a.Re, Im-a.Im  ); }
  __device__ cuComplex operator*(const cuComplex& a) { return cuComplex( Re*a.Re - Im*a.Im, 
                                                                         Im*a.Re + Re*a.Im ); }
  __device__ cuComplex operator*(const float a)      { return cuComplex( Re*a,  Im*a       ); }
};

/*-----------------------------------------------------------------------------------------------*/

__global__ void view_m (cuComplex A[Nc][Nc]);                                 // print matrix

__device__ void copy_m (cuComplex A[Nc][Nc], cuComplex B[Nc][Nc]);       // copy A->B

__device__ void mul_m  (cuComplex A[Nc][Nc],  cuComplex B[Nc][Nc],       // multiply A*B -> C
                                      cuComplex C[Nc][Nc]);
__device__ void add_m  (cuComplex A[Nc][Nc],  cuComplex B[Nc][Nc],       // add      A+B -> C
                                      cuComplex C[Nc][Nc]);
__device__ void sub_m  (cuComplex A[Nc][Nc],  cuComplex B[Nc][Nc],       // subtract A-B -> C
                                      cuComplex C[Nc][Nc]);
__device__ void equ_m  ( cuComplex a,                                         // set A = a*Id
                                      cuComplex A[Nc][Nc]);
__device__ void conj_m                          (cuComplex A[Nc][Nc]);        // A -> conjugate(A)
__device__ void dag_m                           (cuComplex A[Nc][Nc]);        // A -> A^\dagger

__device__ cuComplex det_m                 (cuComplex A[Nc][Nc]);        // returns det(A)

__device__ cuComplex trace                 (cuComplex A[Nc][Nc]);        // returns tr(A)

__device__ void gramschmidt                     (cuComplex A[Nc][Nc]);        // Gram-Schmidt reduce A

__global__ void suN_m                           (cuComplex A[Nc][Nc]);        // project into SU(Nc) 

/*-----------------------------------------------------------------------------------------------*/
