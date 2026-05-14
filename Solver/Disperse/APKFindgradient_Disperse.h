#ifndef __APKFindgradient_Disperse_h
#define __APKFindgradient_Disperse_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"
#include"APKCalcGradVar.h"
#include "APKCalcGradVar_Volfrac.h"


double Findgradient_Disperse(Base *base,Zone *zone,int index[3],int variableTag);

#endif

