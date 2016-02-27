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


void ic(int h) {
  /*
   * h = -1 or +1 to all links
   *      0 gives random alignments
   *
   * general link elements are
   */
  int x[4];
  for (x[0]=0; x[0]<NX; x[0]++)  {
  for (x[1]=0; x[1]<NX; x[1]++)    {
  for (x[2]=0; x[2]<NX; x[2]++)      {
  for (x[3]=0; x[3]<NX; x[3]++)        {

    for (int d=0; d<4; d++) {

                    if (h== 0) link[ x[0] ][ x[1] ][ x[2] ][ x[3] ][d] = ( rand() % 2 )*2-1   ;
                    if (h==+1) link[ x[0] ][ x[1] ][ x[2] ][ x[3] ][d] =  1                   ;
               else if (h==-1) link[ x[0] ][ x[1] ][ x[2] ][ x[3] ][d] = -1                   ;

            }
        }
      }
    }
  }
  return;
}


/*-----------------------------------------------------------------------------------------------*/

void shift( int x[], int d, int steps) {
  /*
   * d = 1...4     (direction)
   * x = steps      (can be neg)
   */
  x[d] = (  x[d]  +  steps +NX ) % NX;
  return;
}

double update(double beta, int x[], int d, int dim) {
  double complex uij;
  double complex S=0., Spl=1.;
  double prob, action; 

  for (int dp=0; dp<dim; dp++)    {
    if (dp!=d)                    {
        // assume this link is +1
        //                              [ redo later for `streamlined code' ]

      shift(x, dp,  -1);
      Spl  = link[ x[0] ][ x[1] ][ x[2] ][ x[3] ][dp];
      Spl *= link[ x[0] ][ x[1] ][ x[2] ][ x[3] ][d ];
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

    }
  } 

  if (zn==0) {
    uij = cexp( I*2*M_PI*((double) (rand()))/( (double) RAND_MAX ));
  }
  else       {
    uij = cexp( I*2*M_PI*((double) (rand()%zn))/(double) zn);
  }


  double dS = -creal( S*(uij - link[ x[0] ][ x[1] ][ x[2] ][ x[3] ][d] ) );
  action = 0.;

  // Boltzmann factor
  prob = exp( - beta*dS );
  /*prob = prob/( prob + (1./prob) );*/

  if ( ( (float) rand() )/RAND_MAX < prob ) {
    link[ x[0] ][ x[1] ][ x[2] ][ x[3] ][d] = uij;
    action += creal( S*uij );
  }
  else {
    action += creal( S*link[ x[0] ][ x[1] ][ x[2] ][ x[3] ][d] );
  }

  return action;
}

double sweep( double beta, int dim ) {

  int       x[4], d;
  double action = 0.0;

  for (int run=0; run<calls; run++)  {

    x[0] = rand() % NX;
    x[1] = rand() % NX;
    x[2] = rand() % NX;
    x[3] = rand() % NX;
    d    = rand() % dim;

    action += update( beta, x, d, dim );
  }

  action /= 2.*( ((double) dim) - 1. )*((double) calls);
  return 1. - action;

}


