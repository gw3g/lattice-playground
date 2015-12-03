/*
 * Author: greg jackson
 * Date: Dec 03 2015
 * Z_N gauge theory on a 4D lattice
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <complex.h>

/*-----------------------------------------------------------------------------------------------*/

/* external parameters */

#define N 10

double complex link[N][N][N][N][4]        ;   // N**4 array
int        calls = 10000000          ;   // MC calls
int        zn    = 0;
           /*srand(time(NULL))    ;*/

/*-----------------------------------------------------------------------------------------------*/

void ic(int h) {
  /*
   * h = -1 or +1 to all links
   *      0 gives random alignments
   *
   * general link elements are
   */
  int x[4];
  for (x[0]=0; x[0]<N; x[0]++)  {
  for (x[1]=0; x[1]<N; x[1]++)    {
  for (x[2]=0; x[2]<N; x[2]++)      {
  for (x[3]=0; x[3]<N; x[3]++)        {

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

void shift( int x[], int d, int steps) {
  /*
   * d = 1...4     (direction)
   * x = steps      (can be neg)
   */
  x[d] = (  x[d]  +  steps +N ) % N;
  /*if (x[d]<0) x[d]+=N;*/
  return;
}

double update(double beta, int x[], int d) {
  double complex uij;
  double S=0., Spl=1.;
  double prob, action; 

  for (int dp=0; dp<4; dp++)    {
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
    uij = cexp( I*2*M_PI*((double) (rand()))/(double) RAND_MAX);
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


  // decide...
  /*if ( ( (float) rand() )/RAND_MAX < prob ) {*/
    /*link[ x[0] ][ x[1] ][ x[2] ][ x[3] ][d ]  =  +1;*/
    /*action += S;*/
    /*action = +( (double) S);*/
  /*}*/
  /*else {*/
    /*link[ x[0] ][ x[1] ][ x[2] ][ x[3] ][d ]  =  -1;*/
    /*action -= S;*/
    /*action = -( (double) S);*/
  /*}*/
  return action;
}

double sweep( double beta ) {

  int       x[4], d;
  double action = 0.0;

  for (int run=0; run<calls; run++)  {

    x[0] = rand() % N;
    x[1] = rand() % N;
    x[2] = rand() % N;
    x[3] = rand() % N;
    d    = rand() % 4;

    action += update( beta, x, d );
  }

  action /= 6.*((double) calls);
  return 1. - action;

}


/*-----------------------------------------------------------------------------------------------*/

FILE *file; char fname[40];

int main() {

  srand(time(NULL))    ;
  ic(1);

  double beta = .0, db = .05, S;

  sprintf(fname, "out/z2, 4D, heat.dat");
  file = fopen(fname, "w+");
  for(beta = 0.0; beta<2.1+db; beta+=db) {
    S = sweep(beta);
    printf(       "%g\t %.8f\n", beta, S );
    fprintf(file, "%g\t %.8f\n", beta, S );
  }
  fclose(file);

  printf("\n\n");

  sprintf(fname, "out/z2, 4D, cool.dat");
  file = fopen(fname, "w+");
  for(beta = 2.1; beta>0.0-db; beta-=db) {
    S = sweep(beta);
    printf(       "%g\t %.8f\n", beta, S );
    fprintf(file, "%g\t %.8f\n", beta, S );
  }
  fclose(file);

  return 0;
}


