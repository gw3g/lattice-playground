#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "matrix.h"

/*-----------------------------------------------------------------------------------------------*/

#define NX 16
#define DIM 4

typedef struct Group { double complex U[Nc][Nc]; } Group;

/*-----------------------------------------------------------------------------------------------*/

extern Group         *ulinks                    ;
extern int            calls                      ;   // MC calls
extern int            zn                         ;   // if 0 --> U(1)

/*-----------------------------------------------------------------------------------------------*/


/*
 *  conventions:      mu =    0, 1, ..., DIM-1     ( direction )
 *                x[DIM] = { t, x, y, z, ... }     ( coordinates )
 */

inline int sites() {
  // --- number of sites ~ NX^DIM
  int o = 1;
  for (int i=0; i<DIM; i++) o *= NX;
  return o;
}

inline int Idx(int x[DIM], int mu) {
  // --- lattice index from coordinates & direction
  int offset=0, nsites=1;
  for (int i=0; i<DIM; i++) {    offset += x[i]*nsites;
                                 nsites *= NX;               };
  return offset + mu*nsites;
}

inline int xMu(int idx, int x[DIM]) {
  // --- returns mu and modifies x[..]
  int tI = idx;
  for (int i=0; i<DIM; i++) {     x[i]  = tI%NX;
                                  tI   /= NX;        };
  return tI;
}

inline void shift_x( int x[DIM], int mu, int steps) {
  // --- shifts coordinates in mu direction
  x[mu] = (x[mu]+steps+NX)%NX;
  return;
}

/*-----------------------------------------------------------------------------------------------*/

// Z_n and U(1) groups ``simple''
Group *init(int);
double update_Zn(double, Group *, int *, int );
double sweep_Zn( double, Group * );

// SU(Nc) group
Group *init_COLD();
Group *init_HOT ();
double update(double, Group *, int *, int );
double sweep(double, Group *);
double Wloop( int,int, Group *);

