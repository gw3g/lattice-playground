#include <stdio.h>
#include <../include/core.h>

/*-----------------------------------------------------------------------------------------------*/

double db = .05;
FILE *file; char fname[40];

void eval_run( double Bi, double Bf, int dim) {
  ic(1);

  double beta, S;

  if (Bi<Bf) {                                                      // COOLING
    sprintf(fname, "out/data/Z%d, %dD, cool.csv", zn, dim);
    file = fopen(fname, "w+");
    for (beta=Bi; beta<Bf+db; beta+=db) 
        {   S = sweep(beta, dim);    
            fprintf(file, "%.8f, %.8f\n", beta, S );
            printf(       "%.8f, %.8f\n", beta, S );    }
  }
  else if (Bi>Bf) {                                                 // HEATING
    sprintf(fname, "out/data/Z%d, %dD, heat.csv", zn, dim);
    file = fopen(fname, "w+");
    for (beta=Bi; beta>Bf-db; beta-=db) 
        {   S = sweep(beta, dim);    
            fprintf(file, "%.8f, %.8f\n", beta, S );
            printf(       "%.8f, %.8f\n", beta, S );    }
  }
  fclose(file);                                                                             return;
}

/*-----------------------------------------------------------------------------------------------*/

int main() {

  srand(time(NULL))    ;

  /*eval_run(.0, 2.1, 4);*/
  /*eval_run(2.1, .0, 4);*/

  eval_run(.0, 2.1, 2);
  eval_run(2.1, .0, 2);

  return 0;
}


