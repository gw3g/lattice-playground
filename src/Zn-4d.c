/*
 * Author: greg jackson
 * Date: Dec 03 2015
 * Z_N gauge theory on a d<=4 lattice
 * lim_{N->inf} is included, U(1)
 *
 */

#include "../include/core.h"
#include "../include/matrix.h"

/*-----------------------------------------------------------------------------------------------*/


Group *init(int h) {
  int nlinks = DIM*sites();
  Group *L = (Group *)malloc(nlinks*sizeof(Group));
  // --- populate with unit matrices
  for (int i=0; i<nlinks; i++) {
         if (h== 0) L[i].U[0][0] = ( rand() % 2 )*2-1   ;
         if (h==+1) L[i].U[0][0] =  1                   ;
    else if (h==-1) L[i].U[0][0] = -1                   ;
  }
  return L;
}

/*-----------------------------------------------------------------------------------------------*/

double update_Zn(double beta, Group *lattice, int x[DIM], int mu ) {
  double complex uij;
  double complex S=0., Spl=1.;
  double prob, action; 
  int l[7]; // list of link indices (e.g. l[1] = globalIdx of l1 )

  for (int nu=0; nu<DIM; nu++)    {
    if (nu!=mu)                    {
        // assume this link is +1
        //                              [ redo later for `streamlined code' ]

    // do a "lap"
                           l[0] = Idx( x, mu );
    shift_x( x, nu, -1);   l[1] = Idx( x, nu ); Spl  =  lattice[l[1]].U[0][0];
                           l[2] = Idx( x, mu ); Spl *=  conj(lattice[l[2]].U[0][0]);
    shift_x( x, mu, +1);   l[3] = Idx( x, nu ); Spl *=  conj(lattice[l[3]].U[0][0]);
    S   += Spl;
    shift_x( x, nu, +1);   l[4] = Idx( x, nu ); Spl  =  conj(lattice[l[4]].U[0][0]);
    shift_x( x, nu, +1);
    shift_x( x, mu, -1);   l[5] = Idx( x, mu ); Spl *=  lattice[l[5]].U[0][0];
    shift_x( x, nu, -1);   l[6] = Idx( x, nu ); Spl *=  lattice[l[6]].U[0][0];
    S   += Spl;
  }
  } 

  if (zn==0) {
    uij = cexp( I*2*M_PI*((double) (rand()))/( (double) RAND_MAX ));
  }
  else       {
    uij = cexp( I*2*M_PI*((double) (rand()%zn))/(double) zn);
  }


  double dS = -creal( S*(uij - lattice[l[0]].U[0][0] ) );
  action = 0.;

  // Boltzmann factor
  prob = exp( - beta*dS );
  /*prob = prob/( prob + (1./prob) );*/

  if ( ( (float) rand() )/RAND_MAX < prob ) {
    lattice[l[0]].U[0][0] = uij;
    action += creal( S*uij );
  }
  else {
    action += creal( S*lattice[l[0]].U[0][0] );
  }

  return action;
}

double sweep_Zn( double beta, Group *lattice ) {

  int       x[DIM], mu, nsites=sites(), nlinks = DIM*nsites;
  double action = 0.0;

  for (int run=0; run<calls; run++) for (int i=0; i<nlinks; i++)
      {    mu = xMu(i, x); action += update_Zn( beta, lattice, x, mu );    }

  action /= 2.*( ((double) DIM) - 1. )*((double) calls)*((double) nlinks);
  return 1. - action;

}


