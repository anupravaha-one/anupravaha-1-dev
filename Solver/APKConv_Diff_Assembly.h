#ifndef __APKConv_Diff_Assembly_h
#define __APKConv_Diff_Assembly_h
#define Greaterof(A,B) (A>=B)?A:B

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"

#include"APKMakeDiffusionCoeffDomain.h"
#include"APKCalcvarCoeffAndPrevCoeff.h"

double Conv_Diff_Assembly(Base *base,Zone *zone,int g,int vartag);

#endif

