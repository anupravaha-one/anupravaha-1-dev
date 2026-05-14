#ifndef __APKCalculateNormalVector_h
#define __APKCalculateNormalVector_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"

int CalculateNormalVector(Base *CgnsBase, Zone *z, int index[3], int face, double *surfVector);
int CalculateVector(Base *CgnsBase, double *CalcX, double *CalcY, double *CalcZ);

#endif
