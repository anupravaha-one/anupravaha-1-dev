#ifndef __APKCalcCoefficients2d_h
#define __APKCalcCoefficients2d_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"

double CalcCoefficients2d(int index[3],double vol_p,double u[18],int *size,double *beta,double *DiaN[18]);

#endif

