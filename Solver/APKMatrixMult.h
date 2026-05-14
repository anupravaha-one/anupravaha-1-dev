#ifndef __APKMatrixMult_h
#define __APKMatrixMult_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"

double MatrixMult(int index[3],int *size,double u_p,double u[18],double *dia,double *DiaN[18],double *LHS,double *RHS);

#endif

