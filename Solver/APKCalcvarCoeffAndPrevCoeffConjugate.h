#ifndef __APKCalcvarCoeffAndPrevCoeffConjugate_h
#define __APKCalcvarCoeffAndPrevCoeffConjugate_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"

double CalcvarCoefficientConjugate(Base *base,Zone *zone,int variableTag,int INDEX);
double CalcvarCoefficientPrevConjugate(Base *base,Zone *zone,int variableTag,int INDEX);

#endif
