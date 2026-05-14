#ifndef __APKMatrixMultCorner_h
#define __APKMatrixMultCorner_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"

double MatrixMultCorner(int index[3],int *size,double u_p,double u[18],double *dia,double *DiaN[18],double *LHS,double *RHS);

#endif

