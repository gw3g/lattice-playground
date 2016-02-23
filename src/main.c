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

double complex link[NX][NX][NX][NX][4]        ;   // NX**4 array
int        calls = 10000                   ;   // MC calls
int        zn    = 2                      ;   // if 0 --> U(1)
matrix     *ulinks                        ;


void eval_Zn  ( double Bi, double Bf, int dim);
void eval_U1  ( double Bi, double Bf, int dim);
void eval_SUn ( double Bi, double Bf);

/*-----------------------------------------------------------------------------------------------*/

int main() {

  srand(time(NULL))    ;

/*
  eval_Zn(.0, 2.1, 4);
  eval_Zn(2.1, .0, 4);

  eval_Zn(.0, 2.1, 3);
  eval_Zn(2.1, .0, 3);

  eval_Zn(.0, 2.1, 2);
  eval_Zn(2.1, .0, 2);
  */
  /*double action;*/
  /*ulinks = init_COLD( &action );*/
  /*ulinks = init_COLD( &action );*/
  /*printf(" %g \n", action);*/

  /*view_m(ulinks[1].U);*/

  eval_SUn(10., .01);
  eval_SUn(.01, 10.);
/*
  double b =5, db=.2, s;
  for (int i=0; i<30; b+=db, i++) {
    s = monte(b, ulinks, &action);
    printf("%g,  %g \n", b, s);
  } */

  /*int x[4] = {1,2,2,2};*/

  /*update(2.,  ulinks, x, 1);*/

  /*calls = 10000000;*/
  /*eval_U1(4., .0, 4);*/
  /*eval_U1(0., 4., 4);*/

  return 0;
}

/*-----------------------------------------------------------------------------------------------*/

/*double db = .05;*/
void therm( double, double ); // --- for thermometer bar
int Nbeta = 20;
FILE *file; char fname[40];

void eval_Zn( double Bi, double Bf, int dim) {
  ic(1);

  double beta=Bi, S, db=(Bf-Bi)/( (double) Nbeta );

  if (Bi<Bf) {                                                      // COOLING
    sprintf(fname, "out/data/Z%d_cool_(d=%d, N=%d).csv", zn, dim, NX);
  }
  else if (Bi>Bf) {                                                 // HEATING
    sprintf(fname, "out/data/Z%d_heat_(d=%d, N=%d).csv", zn, dim, NX);
  }

  file = fopen(fname, "w+");

  fprintf(file,   "# d=%d lattice, w/ group action Z_%d \n", dim, zn                      );
  printf(         "# d=%d lattice, w/ group action Z_%d \n", dim, zn                      );
  fprintf(file,   "# MC calls %d\n", calls                                                );
  fprintf(file,   "#\n"                                                                   );
  fprintf(file,   "# beta,    action  \n"                                                 );

  for (int i=0; i<Nbeta; beta+=db, i++) 
      {   S = sweep(beta, dim);
          fprintf(file, "%.8f, %.8f\n", beta, S );
          printf(       "%.8f, %.8f\n", beta, S );    }
 
  fclose(file);                                                                             return;
}

void eval_U1( double Bi, double Bf, int dim) {
  zn=0;ic(1);

  double beta=Bi, S, db=(Bf-Bi)/( (double) Nbeta );

  if (Bi<Bf) {                                                      // COOLING
    sprintf(fname, "out/data/U(1)_cool_(d=%d, NX=%d).csv", dim, NX);
  }
  else if (Bi>Bf) {                                                 // HEATING
    sprintf(fname, "out/data/U(1)_heat_(d=%d, NX=%d).csv", dim, NX);
  }

  file = fopen(fname, "w+");

  fprintf(file,   "# d=%d lattice, w/ group action U(1) \n", dim                          );
  printf(         "# d=%d lattice, w/ group action U(1) \n", dim                          );
  fprintf(file,   "# MC calls %d\n", calls                                                );
  fprintf(file,   "#\n"                                                                   );
  fprintf(file,   "# beta,    action  \n"                                                 );

  for (int i=0; i<Nbeta; beta+=db, i++) 
      {   S = sweep(beta, dim);
          fprintf(file, "%.8f, %.8f\n", beta, S );
          printf(       "%.8f, %.8f\n", beta, S );    }
 
  fclose(file);                                                                             return;
}

void eval_SUn( double Bi, double Bf) {

  double beta=Bi, S, db=(Bf-Bi)/( (double) Nbeta );
  printf("\n d = %d lattice (NX=%d) w/ SU(%d) gauge group \n", DIM, NX, Nc);

  if (Bi<Bf) {                                                      // COOLING
    sprintf(fname, "out/data/SU(%d)_cool_(d=%d, NX=%d).csv", Nc, DIM, NX);
    ulinks = init_HOT( ); printf("\n :: COOLING :: \n");
  }
  else if (Bi>Bf) {                                                 // HEATING
    sprintf(fname, "out/data/SU(%d)_heat_(d=%d, NX=%d).csv", Nc, DIM, NX);
    ulinks = init_COLD( ); printf("\n :: HEATING :: \n ");
  }

  file = fopen(fname, "w+");

  fprintf(file,   "# d=%d lattice, w/ group action SU(%d) \n",   DIM, Nc                  );
  fprintf(file,   "# MC calls %d\n", calls                                                );
  fprintf(file,   "#\n"                                                                   );
  fprintf(file,   "# beta,    action  \n"                                                 );

  for (int i=0; i<Nbeta; beta+=db, i++) 
      {   S = sweep(beta, ulinks);
          fprintf(file, "%.8f, %.8f\n", beta, S );
          therm(S, beta ); 
          /*printf(       "%.8f, %.8f\n", beta, S );    */
      }; printf("\n");
 
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
