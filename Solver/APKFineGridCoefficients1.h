#ifndef __APKFineGridCoefficients1_h
#define __APKFineGridCoefficients1_h
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

double FineGridCoefficients1(Base *base,Zone *zone,int g,int vartag);
//double Greaterof(double A,double B);

#endif

