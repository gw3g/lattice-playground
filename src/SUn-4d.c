/*
 * Author: greg jackson
 * Date: Feb 18 2016
 * SU(N) gauge theory on a d-dim lattice
 *
 */
#include "../include/core.h"
#include "../include/matrix.h"

/*-----------------------------------------------------------------------------------------------*/
// various Matrix-struct functions defined here

Group rmg() {    // --- random matrix generator
  Group m;
  for (int i=0; i<Nc; i++) {
    for (int j=0; j<Nc; j++) {
      m.U[i][j] = 3.*( rand()/((double) RAND_MAX)-.5 )+3.*I*( rand()/((double) RAND_MAX)-.5 );
    }
    m.U[i][i] += 1.;
  }

  suN_m(m.U);     // --- project onto SU(Nc)
  return m;
}

Group *init_COLD() {
  int nlinks = DIM*sites();
  Group *L = (Group *)malloc(nlinks*sizeof(Group));
  // --- populate with unit matrices
  for (int i=0; i<nlinks; i++)  equ_m(1., L[i].U);
  return L;
}

Group *init_HOT() {
  int nlinks = DIM*sites();
  Group *L = (Group *)malloc(nlinks*sizeof(Group));
  // --- populate with unit matrices
  for (int i=0; i<nlinks; i++)  L[i] = rmg();
  return L;
}

/*-----------------------------------------------------------------------------------------------*/

void    staple( Group *lattice, int x[DIM], int mu, int nu, Group *st );
double  plaq( Group l, Group *st );
double  update( double beta, Group *lattice, int x[DIM], int mu );
double  sweep( double beta, Group *lattice );

/*-----------------------------------------------------------------------------------------------*/

void staple( Group *lattice, int x[DIM], int mu, int nu, // out:
                                                          Group *st) {
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
  Group l1, l2, l3, l4, l5, l6;

  // do a "lap"
                         l[0] = Idx( x, mu );
  shift_x( x, nu, -1);   l[1] = Idx( x, nu ); copy_m( lattice[l[1]].U, l1.U );
                         l[2] = Idx( x, mu ); copy_m( lattice[l[2]].U, l2.U );
  shift_x( x, mu, +1);   l[3] = Idx( x, nu ); copy_m( lattice[l[3]].U, l3.U );
  shift_x( x, nu, +1);   l[4] = Idx( x, nu ); copy_m( lattice[l[4]].U, l4.U );
  shift_x( x, nu, +1);
  shift_x( x, mu, -1);   l[5] = Idx( x, mu ); copy_m( lattice[l[5]].U, l5.U );
  shift_x( x, nu, -1);   l[6] = Idx( x, nu ); copy_m( lattice[l[6]].U, l6.U );

  // --- orientation is CLOCKWISE
  /*dag_m( l2.U ); dag_m( l3.U ); dag_m( l4.U );*/
  dag_m( l2.U ); dag_m( l3.U ); dag_m( l5.U ); dag_m( l6.U );

  // st[0] = bottom                       st[1] = top
     mul_m(    l1.U, l2.U, st[0].U );     mul_m(    l6.U, l5.U, st[1].U );
     mul_m( st[0].U, l3.U, st[0].U );     mul_m( st[1].U, l4.U, st[1].U );

  return;
}

double plaq(Group l, Group *st) {
  Group pl[2]; double S_plaq;

  // --- combine top & bottom staples
  /*for (int z=0; z<2; z++) {mul_m( l.U, st[z].U, pl[z].U ); dag_m( l.U );}*/
  for (int z=0; z<2; z++) {mul_m( l.U, st[z].U, pl[z].U );}

  add_m( pl[0].U, pl[1].U, pl[0].U );
  S_plaq = 2. - creal( trace( pl[0].U ) )/( (double) Nc );

  return S_plaq;
}

/*double S_tot(Group *lattice) {*/
  /*Group st[2], temp[2], l;*/
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

double update(double beta, Group *lattice, int x[DIM], int mu) {
  Group  m   = rmg(),                  // new link variable
         l   = lattice[ Idx(x,mu) ];   // current -- " --
  Group st[2];
  for (int z=0; z<2; z++)  equ_m(.0,st[z].U);

  double S_m = .0, S_l = .0;

  for (int nu=0; nu<DIM; nu++)                            {
    if (nu!=mu)                                           {
                    staple( lattice, x, mu, nu, st );     // loop over plaquette
                                S_m += plaq( m, st );     // add to action
                                S_l += plaq( l, st );     // 
                                                          }
    /*   p_{bot/top}   */                                 }

  double dS = S_m-S_l;
  /*printf( "dS = %.5f \n", S_m );*/

  // Boltzmann factor
  double prob = exp( -beta*dS );
  /*prob = prob/(prob+1./prob);*/

  if ( ( (float) rand() )/RAND_MAX < prob ) { copy_m(m.U, lattice[ Idx(x,mu) ].U ); dS=S_m;}
  else                                      { dS=S_l; }
  return dS;
}

double sweep( double beta, Group *lattice) {

  int       x[DIM], mu, nsites=sites(), nlinks = DIM*nsites;
  double action = 0.0;

  for (int run=0; run<calls; run++) for (int i=0; i<nlinks; i++) 
      {    mu = xMu(i, x); action += update( beta, lattice, x, mu );    }

  action /= 2.*( ((double) DIM) - 1. )*((double) calls*nlinks);
  /*printf(" %g   \n", *action );*/
  return action;
  /*return action;*/
}

double Wloop( int R, int T, Group *lattice ) {
  int mu = 0, nu = 1, nsites=sites();
  Group line;
  int link, x[DIM];
  double complex tr = 0.;
  for (int i=0; i<1; i++) {
    xMu(i, x);
    // consider loop with one corner here
    equ_m(1., line.U);
    for (int r=0; r<R; r++) {
      link = Idx(x,mu); mul_m(line.U, lattice[link].U, line.U); shift_x( x, mu, +1);
    }
    for (int t=0; t<T; t++) {
      link = Idx(x,nu); mul_m(line.U, lattice[link].U, line.U); shift_x( x, nu, +1);
    }
    for (int r=R; r>0; r--) {
      link = Idx(x,mu); mul_m(line.U, lattice[link].U, line.U); shift_x( x, mu, -1);
    }
    for (int t=T; t>0; t--) {
      link = Idx(x,nu); mul_m(line.U, lattice[link].U, line.U); shift_x( x, nu, -1);
    }
    tr += fabs(  creal(trace(line.U))  );
  }
  return creal(tr)/( (double) nsites );
}
