#ifndef __APKRadiationSolutionInt3D_h
#define __APKRadiationSolutionInt3D_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
 
int RadiationSolutionInt3D(Base *base, Zone *zone, double *Intface, double Diag[6], double *NewInt, double *OldInt, int index[5], double alpha);

#endif
