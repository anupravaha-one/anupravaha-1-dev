#ifndef __APKTURBCalcDoubleDerivative2d_h
#define __APKTURBCalcDoubleDerivative2d_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"../Common/APKData_Structures.h"
#include"../Common/APKMakeSlabPtrArray.h"
#include"../Common/APKMakeDomainArray.h"

int DoubleDerivative2d(int index[3],double *domain,double *volumes,int *ptrSize,double *slabVariable[6],double *slabVolume[6],double *beta,double DiffFlux[6]);

#endif
