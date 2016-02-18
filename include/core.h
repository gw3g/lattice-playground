#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <complex.h>

#include "matrix.h"

/*-----------------------------------------------------------------------------------------------*/

#define N 3

typedef struct matrix { double complex U[Nc][Nc]; } matrix;

/*-----------------------------------------------------------------------------------------------*/

extern double complex link[N][N][N][N][4]        ;   // N**4 array
extern matrix         *ulinks                    ;
extern int            calls                      ;   // MC calls
extern int            zn                         ;   // if 0 --> U(1)

/*-----------------------------------------------------------------------------------------------*/

void ic(int);
double update(double, int *, int, int );
double sweep( double, int );

double update(double, matrix *, int *, int );
matrix *init();

double monte(double, matrix *);

