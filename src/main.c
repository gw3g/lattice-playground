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
int        calls = 1000000                ;   // MC calls
int        zn    = 2                      ;   // if 0 --> U(1)

/*-----------------------------------------------------------------------------------------------*/

/*double db = .05;*/
int Nbeta = 20;
FILE *file; char fname[40];

void eval_run( double Bi, double Bf, int dim) {
  ic(1);

  double beta=Bi, S, db=(Bf-Bi)/( (double) Nbeta );

  printf("test\n");

  if (Bi<Bf) {                                                      // COOLING
    sprintf(fname, "out/data/Z%d_cool_(d=%d, N=%d).csv", zn, dim, N);
  }
  else if (Bi>Bf) {                                                 // HEATING
    sprintf(fname, "out/data/Z%d_heat_(d=%d, N=%d).csv", zn, dim, N);
  }

  file = fopen(fname, "w+");
  printf("test\n");

  fprintf(file,   "# d=%d lattice, w/ group action Z_%d \n", dim, zn                      );
  printf(         "# d=%d lattice, w/ group action Z_%d \n", dim, zn                      );
  fprintf(file,   "# MC calls %d\n", calls                                                );
  fprintf(file,   "#\n"                                                                   );
  fprintf(file,   "# beta,    action  \n"                                                 );

  for (int i=0; i<Nbeta; beta+=db, i++) 
      {   S = sweep(beta, dim);
          fprintf(file, "%.8f, %.8f\n", beta, S );
          printf(       "%.8f, %.8f\n", beta, S );    }
 
  /*else if (Bi>Bf) {                                                 // HEATING*/
    /*sprintf(fname, "out/data/Z%d, %dD, N=%d, heat.csv", zn, dim, N);*/
    /*file = fopen(fname, "w+");*/
    /*for (beta=Bi; beta>Bf-db; beta-=db) */
        /*{   S = sweep(beta, dim);    */
            /*fprintf(file, "%.8f, %.8f\n", beta, S );*/
            /*printf(       "%.8f, %.8f\n", beta, S );    }*/
  /*}*/
  fclose(file);                                                                             return;
}

/*-----------------------------------------------------------------------------------------------*/

int main() {

  srand(time(NULL))    ;


  eval_run(.0, 2.1, 4);
  eval_run(2.1, .0, 4);

  eval_run(.0, 2.1, 3);
  eval_run(2.1, .0, 3);

  eval_run(.0, 2.1, 2);
  eval_run(2.1, .0, 2);

  return 0;
}


