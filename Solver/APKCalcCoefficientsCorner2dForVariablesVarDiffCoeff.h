#ifndef __APKCalcCoefficientsCorner2dForVariablesVarDiffCoeff_h
#define __APKCalcCoefficientsCorner2dForVariablesVarDiffCoeff_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeDiffusionCoeffDomain.h"
#include"APKMakeDiffusionCoeffSlabPtrArray.h"
#include"APKCalcCoefficientsCorner2dForVariablesVarDiffCoeffMPFV.h"

double CalcCoefficientsCorner2dForVariablesVarDiffCoeff(Zone *zone,int index[3],double vol_p,double u[18],int *size,double *beta,double *DiaN[18],int vartag);

#endif

