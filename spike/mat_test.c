#include "../include/matrix.h"
#include <stdlib.h>
#include <time.h>

#define nl printf("\n")


int main() {
  double complex xx[Nc][Nc], yy[Nc][Nc], zz[Nc][Nc] = { .0 };

  nl; printf("  equ_m( 2.5, yy ) = "); nl; nl;
  equ_m(2.5,yy); view_m(yy);
  printf("  det yy = %g + i %g", creal(det_m(yy)), cimag(det_m(yy)) ); nl;
  printf("  tr yy = %g + i %g", creal(trace(yy)), cimag(trace(yy)) ); nl;

  // define matrix xx
  for (int i=0; i<Nc; i++)
    for (int j=0; j<Nc; j++)
      xx[i][j] = 1.+2.*i + j*(j+.5*i)*I;
  add_m(yy,xx,xx);

  nl; printf(" xx = "); nl; nl;
  view_m(xx); nl;

  // get adventurous with zz ...
  printf("Adding matrices: "); nl;
  add_m( xx, yy, zz );
  view_m(xx);
  printf("+");
  view_m(yy);
  printf("=");
  view_m(zz); nl;

  // check some determinants:
  printf("Some determinants :"); nl; nl;
  printf("  det xx = %g + i %g \n", creal(det_m(xx)), cimag(det_m(xx)) );
  printf("  det yy = %g + i %g \n", creal(det_m(yy)), cimag(det_m(yy)) ); 
  printf("  det zz = %g + i %g \n", creal(det_m(zz)), cimag(det_m(zz)) ); nl;

  printf("Then we can turn them into SU(Nn) matrices: "); nl; nl;
  suN_m(xx);
  suN_m(yy);
  suN_m(zz);
  view_m(zz);
  view_m(yy);
  view_m(zz);

  printf("det xx = %g + i %g \n", creal(det_m(xx)), cimag(det_m(xx)) );
  printf("det yy = %g + i %g \n", creal(det_m(yy)), cimag(det_m(yy)) );
  printf("det zz = %g + i %g \n", creal(det_m(zz)), cimag(det_m(zz)) ); nl;

  FILE *file; char fname[40];
  sprintf(fname, "out/data/tr_U, Nc=%d.csv", Nc);
  file = fopen(fname, "w+");
  srand(time(NULL));
  fprintf(file, "# trace of 'randomly' generated SU(%d) matrices\n", Nc);
  fprintf(file, "#\n");

  double complex aa[Nc][Nc], tr;
  for (int l=1; l<2000; l++) {
    for (int i=0; i<Nc; i++)
      for (int j=0; j<Nc; j++)
        aa[i][j] = ( rand()/((double) RAND_MAX)-.5 ) + I*( rand()/((double) RAND_MAX)-.5 );

    suN_m(aa);
    tr = trace(aa);
    /*printf( "%.8f, %.8f\n", creal(tr), cimag(tr) );*/
    fprintf( file, "%.8f, %.8f\n", creal(tr), cimag(tr) );
  }

  return 0;
}
