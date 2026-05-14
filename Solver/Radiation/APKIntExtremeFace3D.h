#ifndef __APKIntExtremeFace3D_h
#define __APKIntExtremeFace3D_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"

double IntExtremeFace3D(Base *base, Zone *zone, double *Intface, double *NewInt, double *OldInt, double *slabInt[6], double *slabvolume[6], int index[5], double alpha, int face);

#endif
