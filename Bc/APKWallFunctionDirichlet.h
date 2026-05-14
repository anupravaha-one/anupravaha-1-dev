#ifndef __APKWallFunctionDirichlet_h
#define __APKWallFunctionDirichlet_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"

double WallFunctionDirichlet(Zone *z,ZoneBoundary *zBc,int vartag,int slabcellindex,double Bv,int flagMPFV);

#endif

