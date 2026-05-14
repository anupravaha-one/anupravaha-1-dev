#ifndef __APKTURBCalcGradSquare_h
#define __APKTURBCalcGradSquare_h

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"
#include"APKCalcGradVar.h"

double CalcGradSquare(Base *base,Zone *zone,int index[3],int variableTag);

#endif

