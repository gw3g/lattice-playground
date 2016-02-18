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

double update(double beta, matrix *lattice, int x[DIM], int mu) {
  matrix m = rnd();

  double complex stpl[Nc][Nc];
  double complex p_t[Nc][Nc];
  double complex p_b[Nc][Nc];
  equ_m(0., p_t);
  equ_m(0., p_b);

  matrix ucurr;

  for (int nu=0; nu<DIM; nu++) {
    if (nu!=mu)                {
      /*
       *    x---->---x
       *    |        |
       *    ^        v
       *    |  <     |
       *    0---?----x
       *    |    >   |
       *    ^        v      ^ nu
       *    |        |      |
       *    x---<----x       --- > mu
       *
       */

      // bottom staple:
      equ_m(0., stpl);
      move(x, nu,  -1);

      ucurr = lattice[ Idx(x,nu) ];
      add_m( stpl, ucurr.U, stpl);

      ucurr = lattice[ Idx(x,mu) ];
      dag_m(ucurr.U);

      mul_m( stpl, ucurr.U, stpl);
      dag_m(ucurr.U);

      move(x, mu, +1);
      ucurr = lattice[ Idx(x,nu) ];
      dag_m(ucurr.U);
      mul_m( stpl, ucurr.U, stpl);
      dag_m(ucurr.U);

      add_m( stpl, p_b, p_b );

      // top staple:
      equ_m(0., stpl);
      move(x, nu,  +1);

      ucurr = lattice[ Idx(x,nu) ];

      printf("ucurr = \n");
      view_m( ucurr.U );

      dag_m(ucurr.U);
      add_m( stpl, ucurr.U, stpl);
      dag_m(ucurr.U);
      printf("ucurr = \n");
      view_m( ucurr.U );

      printf("stpl = \n");
      view_m( stpl );

      move(x, nu,  +1);
      move(x, mu,  -1);
      ucurr = lattice[ Idx(x,mu) ];
      mul_m( stpl, ucurr.U, stpl);

      move(x, nu, -1);
      ucurr = lattice[ Idx(x,nu) ];
      mul_m( stpl, ucurr.U, stpl);

      add_m( stpl, p_t, p_t );

  }
  }

  printf("ptop = \n");
  view_m( p_t );
  printf("pbot = \n");
  view_m( p_b );

  sub_m( lattice[ Idx(x,mu) ].U, m.U, stpl ); // treat stpl as temp mat
  mul_m(p_b, stpl, p_b);
  conj_m(stpl);
  mul_m(p_t, stpl, p_t);
  add_m(p_t, p_b, stpl);
  view_m( stpl );

  double dS = -creal( trace( stpl ) ) / 3.;
  /*printf( "%.5f \n", dS );*/

  // Boltzmann factor
  double prob = exp( - beta*dS );
  /*printf( "%.5f \n", prob );*/
  /*prob = prob/( prob + (1./prob) );*/

  if ( ( (float) rand() )/RAND_MAX < prob ) {
    copy_m(m.U, lattice[ Idx(x,mu) ].U );
    /*view_m(lattice[ Idx(x,mu) ].U);*/
  }
  /*else {*/
    /*action += creal( S*link[ x[0] ][ x[1] ][ x[2] ][ x[3] ][d] );*/
  /*}*/

  printf(" %g   \n", dS );
  /*return action;*/

  return dS;
}

double monte( double beta, matrix *lattice ) {

  int       x[4], mu;
  double action = 0.0;

  for (int run=0; run<calls; run++)  {

    x[0] = rand() % N;
    x[1] = rand() % N;
    x[2] = rand() % N;
    x[3] = rand() % N;
    mu   = rand() % DIM;

    action += update( beta, lattice, x, mu );
    /*printf(" %g   \n", action );*/
  }

  action /= 2.*( ((double) DIM) - 1. )*((double) calls);
  printf(" %g   \n", action );
  return 1. - action;

}


