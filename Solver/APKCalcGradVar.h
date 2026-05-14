#ifndef __APKCalcGradVar_h
#define __APKCalcGradVar_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"
#include"APKMakeDomainArray.h"

double CalcGradVar(Base *base,Zone *zone,double *variable,double *slabVariable[6],int INDEX[3],double D[3]);
double CalcGradViscoeff(Base *base,Zone *zone,double *variable,double *slabVariable[6],int INDEX[3],double D[3]);

#endif

