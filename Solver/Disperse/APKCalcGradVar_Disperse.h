#ifndef __APKCalcGradVar_Disperse_h
#define __APKCalcGradVar_Disperse_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"
#include"APKMakeDomainArray.h"

double CalcGradVar_Disperse(Base *base,Zone *zone,double *variable,double *slabVariable[6],int INDEX[3],double D[3],double *U_g,double *SlabU_g[6]);

#endif

