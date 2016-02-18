/*
 * Author: greg jackson
 * Date: Feb 18 2016
 * SU(N) gauge theory on a d-dim lattice
 *
 */
#include "../include/core.h"
#include "../include/matrix.h"

#define DIM 4

/*-----------------------------------------------------------------------------------------------*/

inline int sites() {
  int o = 1;
  for (int i=0; i<DIM; i++) o *= N;
  return o;
}

int Idx(int x[DIM], int mu) {
  int offset=0, nsites=1;
  for (int i=0; i<DIM; i++) {
    offset += x[i]*nsites;
    nsites *= N;
  }
  return offset + mu*nsites;
}

matrix *init() {
  /*int x[DIM];*/
  int nsites = sites();

  int nlinks = DIM*nsites;
  /*int nplaquettes = DIM*(DIM-1)*nsites/2;*/
  matrix *L = (matrix *)malloc(nlinks*sizeof(matrix));

  for (int i=0; i<nlinks; i++) {
    // fill the lattice with unit matrices
    equ_m(1., L[i].U);
  }
  return L;
}

matrix rnd() {
  srand(time(NULL));
  matrix m;
  for (int i=0; i<Nc; i++)
    for (int j=0; j<Nc; j++)
      m.U[i][j] = ( rand()/((double) RAND_MAX)-.5 ) + I*( rand()/((double) RAND_MAX)-.5 );
  suN_m(m.U);
  return m;
}


void move( int x[], int d, int steps) {
  /*
   * d = 1...4     (direction)
   * x = steps      (can be neg)
   */
  x[d] = (  x[d]  +  steps +N ) % N;
  return;
}

double update(double beta, int x[DIM], int mu) {
  matrix m = rnd();
  complex double staple[Nc][Nc];

  for (int nu=0; nu<DIM; nu++) {
    if (nu!=mu)                {

      int t=1; /*
      shift(x, dp,  -1);
      staple  = ulinks[Idx(x,nu)].U;
      staple  = mul( ulinks[ Idx(, ulinks[ Idx(x,mu
      stpl *= ulinkx[Idx(x,d )];
      shift(x, d,  +1);
      Spl *= link[ x[0] ][ x[1] ][ x[2] ][ x[3] ][dp];

      S   += Spl;

      shift(x, dp, +1);
      Spl  = link[ x[0] ][ x[1] ][ x[2] ][ x[3] ][dp];
      shift(x, dp, +1);
      shift(x, d,  -1);
      Spl *= link[ x[0] ][ x[1] ][ x[2] ][ x[3] ][d ];
      shift(x, dp, -1);
      Spl *= link[ x[0] ][ x[1] ][ x[2] ][ x[3] ][dp];

      S   += Spl;
      */
  }
  }
 
  return .1;
}

