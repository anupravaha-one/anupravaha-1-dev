#ifndef __APKCalcGradVar_Volfrac_h
#define __APKCalcGradVar_Volfrac_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"
#include"APKMakeDomainArray.h"

double CalcGradVar_Volfrac(Base *base,Zone *zone,double *variable,double *slabVariable[6],int INDEX[3],double D[3],int variableTag);

#endif

