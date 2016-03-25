/*
 * Author: greg jackson
 * Date: Dec 02 2015
 * simple gauge theories on d = {2,3,4} lattice
 *
 *
 */

#include <stdio.h>
#include <../include/core.h>

/*-----------------------------------------------------------------------------------------------*/

/* external parameters */

int        calls = 10                     ;   // MC calls
int        zn    = 2                       ;   // if 0 --> U(1)
Group     *ulinks                          ;   // the lattice


void eval_Zn  ( double Bi, double Bf);
void eval_U1  ( double Bi, double Bf);
void eval_SUn ( double Bi, double Bf);

void iter_SUn( double beta, int cls);
void wilson(   double beta );

/*-----------------------------------------------------------------------------------------------*/

int main() {

  srand(time(NULL))    ;
  /*iter_SUn( .2, 10);*/
  wilson(.2);
  /*
  eval_Zn(.0, 1.);
  eval_Zn(1., .0);

  zn = 3;
  eval_Zn(.0, 2.);
  eval_Zn(2., .0);
  zn = 4;
  eval_Zn(.0, 2.);
  eval_Zn(2., .0);

  zn = 5;
  eval_Zn(.0, 3.);
  eval_Zn(3., .0);
  zn = 6;
  eval_Zn(.0, 3.);
  eval_Zn(3., .0);

  zn = 7;
  eval_Zn(.0, 4.);
  eval_Zn(4., .0);
  zn = 8;
  eval_Zn(.0, 4.);
  eval_Zn(4., .0);

  eval_U1(.0, 4.);
  eval_U1(4., .0);
  */

  /*eval_SUn(10., .0);*/
  /*eval_SUn(.0, 10.);*/
/*
  double b =5, db=.2, s;
  for (int i=0; i<30; b+=db, i++) {
    s = monte(b, ulinks, &action);
    printf("%g,  %g \n", b, s);
  } */

  return 0;
}

/*-----------------------------------------------------------------------------------------------*/

/*double db = .05;*/
void therm( double, double ); // --- for thermometer bar
int Nbeta = 40;
FILE *file; char fname[40];

void eval_Zn( double Bi, double Bf) {
  ulinks = init(1);

  double beta=Bi, S, db=(Bf-Bi)/( (double) Nbeta );
  printf("\nDETAILS: d = %d lattice (NX=%d) w/ Z_%d gauge group \n", DIM, NX, zn);

  if (Bi<Bf) {                                                      // COOLING
    sprintf(fname, "out/data/Z%d_cool_(d=%d, NX=%d).csv", zn, DIM, NX);
    printf("\n :: COOLING :: \n\n");
  }
  else if (Bi>Bf) {                                                 // HEATING
    sprintf(fname, "out/data/Z%d_heat_(d=%d, NX=%d).csv", zn, DIM, NX);
    printf("\n :: HEATING :: \n\n");
  }

  file = fopen(fname, "w+");

  fprintf(file,   "# d=%d lattice, w/ group action Z_%d \n", DIM, zn                      );
  fprintf(file,   "# MC calls %d\n", calls                                                );
  fprintf(file,   "#\n"                                                                   );
  fprintf(file,   "# beta,    action  \n"                                                 );

  for (int i=0; i<Nbeta+1; i++) 
      {   S = sweep_Zn(beta, ulinks);
          therm(S, beta ); 
          fprintf(file, "%.8f, %.8f\n", beta, S ); beta+=db;
          /*printf(       "%.8f, %.8f\n", beta, S );    }*/
      }
 
  fclose(file);                                                                             return;
}

void eval_U1( double Bi, double Bf) {
  zn=0; ulinks = init(1);

  double beta=Bi, S, db=(Bf-Bi)/( (double) Nbeta );
  printf("\nDETAILS: d = %d lattice (NX=%d) w/ U(1) gauge group \n", DIM, NX);

  if (Bi<Bf) {                                                      // COOLING
    sprintf(fname, "out/data/U(1)_cool_(d=%d, NX=%d).csv", DIM, NX);
    printf("\n :: COOLING :: \n\n");
  }
  else if (Bi>Bf) {                                                 // HEATING
    sprintf(fname, "out/data/U(1)_heat_(d=%d, NX=%d).csv", DIM, NX);
    printf("\n :: HEATING :: \n\n");
  }

  file = fopen(fname, "w+");

  fprintf(file,   "# d=%d lattice, w/ group action U(1) \n", DIM                          );
  fprintf(file,   "# MC calls %d\n", calls                                                );
  fprintf(file,   "#\n"                                                                   );
  fprintf(file,   "# beta,    action  \n"                                                 );

  for (int i=0; i<Nbeta+1; i++) 
      {   S = sweep_Zn(beta, ulinks);
          therm(S, beta ); 
          fprintf(file, "%.8f, %.8f\n", beta, S ); beta+=db;
          /*printf(       "%.8f, %.8f\n", beta, S );    }*/
      }
 
  fclose(file);                                                                             return;
}

void eval_SUn( double Bi, double Bf) {

  double beta=Bi, S, db=(Bf-Bi)/( (double) Nbeta );
  printf("\nDETAILS: d = %d lattice (NX=%d) w/ SU(%d) gauge group \n", DIM, NX, Nc);

  if (Bi<Bf) {                                                      // COOLING
    sprintf(fname, "out/data/SU(%d)_cool_(d=%d, NX=%d).csv", Nc, DIM, NX);
    ulinks = init_HOT( ); printf("\n :: COOLING :: \n\n");
  }
  else if (Bi>Bf) {                                                 // HEATING
    sprintf(fname, "out/data/SU(%d)_heat_(d=%d, NX=%d).csv", Nc, DIM, NX);
    ulinks = init_HOT( ); printf("\n :: HEATING :: \n\n");
  }

  file = fopen(fname, "w+");

  fprintf(file,   "# d=%d lattice, w/ group action SU(%d) \n",   DIM, Nc                  );
  fprintf(file,   "# MC calls %d\n", calls                                                );
  fprintf(file,   "#\n"                                                                   );
  fprintf(file,   "# beta,    action  \n"                                                 );

  for (int i=0; i<Nbeta+1; i++) 
      {   S = sweep(beta, ulinks);
          fprintf(file, "%.8f, %.8f\n", beta, S );
          therm(S, beta ); beta+=db;
          /*printf(       "%.8f, %.8f\n", beta, S );    */
      }; printf("\n");
 
  fclose(file); free(ulinks);                                                               return;
}

void iter_SUn( double beta, int cls) {

  printf("\nDETAILS: d = %d lattice (NX=%d) w/ SU(%d) gauge group \n", DIM, NX, Nc);

  sprintf(fname, "out/data/SU(%d)_beta=%.2f_(d=%d, NX=%d).csv", Nc, beta, DIM, NX);
  ulinks = init_COLD( );

  file = fopen(fname, "w+");

  fprintf(file,   "# d=%d lattice, w/ group action SU(%d) \n",   DIM, Nc                  );
  fprintf(file,   "#\n"                                                                   );
  fprintf(file,   "# iter,    action  \n"                                                 );

  double S;
  for (int i=0; i<cls; i++) 
      {   S = sweep(beta, ulinks);
          fprintf(file, "%d, %.8f\n", i, S );
          therm(S, beta ); 
          /*printf(       "%.8f, %.8f\n", beta, S );    */
      }; printf("\n");

  printf( "W(3,4) = %g\n", Wloop(3,4,ulinks) );
 
  fclose(file); free(ulinks);                                                               return;
}

void wilson( double beta ) {

  printf("\nDETAILS: d = %d lattice (NX=%d) w/ SU(%d) gauge group \n", DIM, NX, Nc);

  sprintf(fname, "out/data/WILSON_beta=%.2f_(d=%d, NX=%d).csv", beta, DIM, NX);
  ulinks = init_HOT( );

  file = fopen(fname, "w+");

  fprintf(file,   "# d=%d lattice, w/ group action SU(%d) \n",   DIM, Nc                  );
  fprintf(file,   "#\n"                                                                   );
  fprintf(file,   "# R, T, tr(W)      \n"                                                 );

  double S;
  for (int i=0; i<20; i++) 
      {   S = sweep(beta, ulinks);
          therm(S, beta ); 
      };  printf("\n");

  double w;
  for (int R=1; R<NX; R++) for (int T=1; T<NX; T++) {
    w = Wloop( R, T, ulinks );
    printf( "W(%d,%d) = %g\n",      R, T, w );
    fprintf(file, "%d, %d, %.8f\n", R, T, w );
  }
 
  fclose(file); free(ulinks);                                                               return;
}



/*-----------------------------------------------------------------------------------------------*/

void therm( double e, double b ) {
  int width = 60;
  printf(" beta = %.1f : [", b);
  for (int i=0; i<width;++i) {
    int pos = e*width / 1.2;
         if (i<pos)  printf("-");
    else if (i==pos) printf("O");
    else             printf(" ");
  }
  printf("]  %.4f\r", e);
  fflush(stdout);
}
