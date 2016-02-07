#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <complex.h>

/*-----------------------------------------------------------------------------------------------*/

#define N 10

/*-----------------------------------------------------------------------------------------------*/

extern double complex link[N][N][N][N][4]        ;   // N**4 array
extern int            calls                      ;   // MC calls
extern int            zn                         ;   // if 0 --> U(1)

/*-----------------------------------------------------------------------------------------------*/

void ic(int);
double update(double, int *, int, int );
double sweep( double, int );
