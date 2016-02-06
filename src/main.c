#include <stdio.h>
#include <../include/core.h>

/*-----------------------------------------------------------------------------------------------*/

double db = .05;
FILE *file; char fname[40];

void eval_cool( double Bi, double Bf) {
  ic(1);

  double beta = .0, db = .05, S;

  sprintf(fname, "out/data/Z%d, 4D, cool.dat", zn);
  file = fopen(fname, "w+");
  for(beta = 0.0; beta<2.1+db; beta+=db) {
    S = sweep(beta);
    printf(       "%g\t %.8f\n", beta, S );
    fprintf(file, "%g\t %.8f\n", beta, S );
  }
  fclose(file);                                                                             return;
}

/*-----------------------------------------------------------------------------------------------*/

int main() {

  srand(time(NULL))    ;
  ic(1);

  double beta = .0, db = .05, S;

  sprintf(fname, "out/data/z%d, 4D, heat.dat", zn);
  file = fopen(fname, "w+");
  for(beta = 0.0; beta<2.1+db; beta+=db) {
    S = sweep(beta);
    printf(       "%g\t %.8f\n", beta, S );
    fprintf(file, "%g\t %.8f\n", beta, S );
  }
  fclose(file);

  printf("\n\n");

  sprintf(fname, "out/data/z%d, 4D, cool.dat", zn);
  file = fopen(fname, "w+");
  for(beta = 2.1; beta>0.0-db; beta-=db) {
    S = sweep(beta);
    printf(       "%g\t %.8f\n", beta, S );
    fprintf(file, "%g\t %.8f\n", beta, S );
  }
  fclose(file);

  return 0;
}


