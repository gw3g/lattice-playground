#include <stdio.h>
#include <math.h>
#include <complex.h>

/*-----------------------------------------------------------------------------------------------*/

#define Nc 4

/*-----------------------------------------------------------------------------------------------*/

void view_m(double complex  X[Nc][Nc]);                               // print matrix

void mul_m(double complex A[Nc][Nc],  double complex B[Nc][Nc],       // multiply A*B -> C
                                      double complex C[Nc][Nc]);
void add_m(double complex A[Nc][Nc],  double complex B[Nc][Nc],       // add      A+B -> C
                                      double complex C[Nc][Nc]);
void sub_m(double complex A[Nc][Nc],  double complex B[Nc][Nc],       // subtract A-B -> C
                                      double complex C[Nc][Nc]);
void equ_m( double complex a,                                         // set A = a*Id
                                      double complex A[Nc][Nc]);
void conj_m(double complex A[Nc][Nc]);                                // A -> conjugate(A)
void dag_m(double complex A[Nc][Nc]);                                 // A -> A^\dagger

double complex det_m(double complex A[Nc][Nc]);                       // returns det(A)

double complex trace(double complex A[Nc][Nc]);                       // returns tr(A)

void gramschmidt(double complex A[Nc][Nc]);                           // Gram-Schmidt reduce A

void suN_m(double complex A[Nc][Nc]);                                 // project into SU(Nc) 

