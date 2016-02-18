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

double complex link[N][N][N][N][4]        ;   // N**4 array
int        calls = 1                      ;   // MC calls
int        zn    = 2                      ;   // if 0 --> U(1)
matrix     *ulinks                        ;

/*-----------------------------------------------------------------------------------------------*/

/*double db = .05;*/
int Nbeta = 100;
FILE *file; char fname[40];

void eval_Zn( double Bi, double Bf, int dim) {
  ic(1);

  double beta=Bi, S, db=(Bf-Bi)/( (double) Nbeta );

  if (Bi<Bf) {                                                      // COOLING
    sprintf(fname, "out/data/Z%d_cool_(d=%d, N=%d).csv", zn, dim, N);
  }
  else if (Bi>Bf) {                                                 // HEATING
    sprintf(fname, "out/data/Z%d_heat_(d=%d, N=%d).csv", zn, dim, N);
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
    sprintf(fname, "out/data/U(1)_cool_(d=%d, N=%d).csv", dim, N);
  }
  else if (Bi>Bf) {                                                 // HEATING
    sprintf(fname, "out/data/U(1)_heat_(d=%d, N=%d).csv", dim, N);
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
  ulinks = init();

  view_m(ulinks[1].U);

  monte(10.1, ulinks);

  /*int x[4] = {1,2,2,2};*/

  /*update(2.,  ulinks, x, 1);*/

  /*calls = 10000000;*/
  /*eval_U1(4., .0, 4);*/
  /*eval_U1(0., 4., 4);*/

  return 0;
}


