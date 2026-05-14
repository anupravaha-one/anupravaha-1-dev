#ifndef __APKCalcvarCoeffAndPrevCoeff_h
#define __APKCalcvarCoeffAndPrevCoeff_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"

double CalcvarCoefficient(Base *base,Zone *zone,int variableTag,int INDEX);
double CalcvarCoefficientPrev(Base *base,Zone *zone,int variableTag,int INDEX);


#endif

