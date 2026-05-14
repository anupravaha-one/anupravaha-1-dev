#ifndef __APKCalculateDmF_h
#define __APKCalculateDmF_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"

int CalculateDmF(Base *CgnsBase, Zone *z, int index[3], int face, double *surfVector);
int CalculateOmega(Base *CgnsBase, double *CalcOmega);

#endif
