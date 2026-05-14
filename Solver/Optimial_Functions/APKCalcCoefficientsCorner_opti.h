#ifndef __APKCalcCoefficientsCorner_opti_h
#define __APKCalcCoefficientsCorner_opti_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"

double CalcCoefficientsCorner_opti(int index[3],double vol_p,double u[18],int *size,double *beta,double *DiaN[18],double *D[36],double *DN[24],double *Diag[6]);

#endif

