#ifndef __APKIntExtremeCell3D_h
#define __APKIntExtremeCell3D_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
 
double IntExtremeCell3D(Base *base, Zone *zone, double *Intface, double *NewInt, double *OldInt, int index[5], double alpha, int face);

#endif
