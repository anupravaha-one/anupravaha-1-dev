#ifndef __APKWallFunctionDirichletConjugate_h
#define __APKWallFunctionDirichletConjugate_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"

double WallFunctionDirichletConjugate(Zone *z,Interface *IfacePtr,int Face,int vartag,int slabcellindex,int domaincellindex,double Bv);

#endif

