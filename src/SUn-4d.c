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

inline int Idx(int x[DIM], int mu) {
  int offset=0, nsites=1;
  for (int i=0; i<DIM; i++) {
    offset += x[i]*nsites;
    nsites *= N;
  }
  return offset + mu*nsites;
}

void shift_x( int x[DIM], int mu, int steps) {
  /*
   * mu = 1...4     (direction)
   */
  x[mu] = (x[mu]+steps+N)%N;
  return;
}


/*-----------------------------------------------------------------------------------------------*/
// various Matrix-struct functions defined here

matrix rmg() {
  // random matrix generator
  matrix m;
  for (int i=0; i<Nc; i++)  for (int j=0; j<Nc; j++)
      m.U[i][j] = ( rand()/((double) RAND_MAX)-.5 ) + I*( rand()/((double) RAND_MAX)-.5 );
  suN_m(m.U); // project to SU(Nc)
  return m;
}

matrix *init_COLD(double *action) {
  /*int x[DIM];*/
  int nsites = sites();

  int nlinks = DIM*nsites;
  /*int nplaquettes = DIM*(DIM-1)*nsites/2;*/
  matrix *L = (matrix *)malloc(nlinks*sizeof(matrix));

  for (int i=0; i<nlinks; i++) {
    // fill the lattice with unit matrices
    equ_m(1., L[i].U);
  }
  *action = 0.0;
  return L;
}

matrix *init_HOT(double *action) {
  /*int x[DIM];*/
  int nsites = sites();

  int nlinks = DIM*nsites;
  /*int nplaquettes = DIM*(DIM-1)*nsites/2;*/
  matrix *L = (matrix *)malloc(nlinks*sizeof(matrix));

  for (int i=0; i<nlinks; i++) {
    // fill the lattice with unit matrices
    L[i] = rmg();
  }
  *action = 0.0;
  return L;
}

/*-----------------------------------------------------------------------------------------------*/
// MC sampling

void staple( matrix *lattice, int x[DIM], int mu, int nu,
                                            matrix *st) {
      /*
       *              l5            calculate product of two staples
       *          5---->---4        adjoining l0 (not including l0).
       *          |        |        store result in st[2]. Two direc-
       *       l6 ^        ^ l4     tions mu and nu given... looped
       *          |   l0   |        over in calc for Wilson action.
       *          0--->----3
       *          |        |
       *       l1 ^        ^ l3     ^ nu
       *          |   l2   |        |
       *          1--->----2        --- > mu
       */

  int l[7]; // list of link indices (e.g. l[1] = globalIdx of l1 )
  matrix l1, l2, l3, l4, l5, l6;

  // do a "lap"
                         l[0] = Idx( x, mu );
  shift_x( x, nu, -1);   l[1] = Idx( x, nu ); copy_m( lattice[l[1]].U, l1.U );
                         l[2] = Idx( x, mu ); copy_m( lattice[l[2]].U, l2.U );
  shift_x( x, mu, +1);   l[3] = Idx( x, nu ); copy_m( lattice[l[3]].U, l3.U );
  shift_x( x, nu, +1);   l[4] = Idx( x, nu ); copy_m( lattice[l[4]].U, l4.U );
  shift_x( x, nu, +1);
  shift_x( x, mu, -1);   l[5] = Idx( x, mu ); copy_m( lattice[l[5]].U, l5.U );
  shift_x( x, nu, -1);   l[6] = Idx( x, nu ); copy_m( lattice[l[6]].U, l6.U );

  // orientation is CLOCKWISE
  dag_m( l2.U ); dag_m( l3.U ); dag_m( l4.U );

  // st[0] = bottom                       st[1] = top
     mul_m(    l1.U, l2.U, st[0].U );     mul_m(    l6.U, l5.U, st[1].U );
     mul_m( st[0].U, l3.U, st[0].U );     mul_m( st[1].U, l4.U, st[1].U );

  return;
}

double plaq(matrix l, matrix *st) {
  matrix pl[2]; double S_plaq;

  for (int z=0; z<2; z++) {mul_m( l.U, st[z].U, pl[z].U ); dag_m( l.U );}

  add_m( pl[0].U, pl[1].U, pl[0].U );
  S_plaq = creal( trace( pl[0].U ) )/( (double) Nc );

  return S_plaq;
}

/*double S_tot(matrix *lattice) {*/
  /*matrix st[2], temp[2], l;*/
    /*l = lattice[ Idx(x,mu) ];*/
  /*for (int z=0; z<2; z++)  equ_m(0,st[z].U);*/
  /*double neigh = DIM*(DIM-1)/2.;*/

  /*for (int mu=0; mu<DIM; mu++)                            {  l = lattice[ Idx(x,mu) ];*/
  /*for (int nu=0; nu<DIM; nu++)                            {*/
    /*if (nu!=mu)                                           {*/
                    /*staple( lattice, x, mu, nu, temp );   // loop over plaquette*/
                  /*add_m( temp[0].U, st[0].U, st[0].U );   // contributions to*/
                  /*add_m( temp[1].U, st[1].U, st[1].U );   // action ...*/
                               /*S_l = plaq(l, st)/neigh;   // printf(" S_l = %g \n", S_l);*/
                                                          /*}*/
    /*[> p_{bot/top} <]                                     }*/
                                                          /*}*/


/*}*/

double update(double beta, matrix *lattice, int x[DIM], int mu) {
  matrix m = rmg(); matrix l = lattice[ Idx(x,mu) ];
  matrix st[2], temp[2];
  for (int z=0; z<2; z++)  equ_m(0,st[z].U);

  double S_m, S_l;

  for (int nu=0; nu<DIM; nu++)                            {
    if (nu!=mu)                                           {
                    staple( lattice, x, mu, nu, temp );   // loop over plaquette
                  add_m( temp[0].U, st[0].U, st[0].U );   // contributions to
                  add_m( temp[1].U, st[1].U, st[1].U );   // action ...
                                                          }
    /* p_{bot/top} */                                     }

  /*double neigh = (DIM-1);*/
  S_m = 1. - plaq(m, st);                           // printf(" S_m = %g \n", S_m);
  S_l = 1. - plaq(l, st);                           // printf(" S_l = %g \n", S_l);

  double dS = S_l-S_m;
  /*printf( "dS = %.5f \n", S_m );*/

  // Boltzmann factor
  double prob = exp( - beta*dS );

  if ( ( (float) rand() )/RAND_MAX < prob ) { copy_m(m.U, lattice[ Idx(x,mu) ].U ); dS=S_m;}
  else                                      { dS=S_l; }
  return dS;
}

double monte( double beta, matrix *lattice) {

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
  /*printf(" %g   \n", *action );*/
  return 1. - action;
}

