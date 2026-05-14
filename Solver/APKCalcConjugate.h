#ifndef __APKCalcConjugate_h
#define __APKCalcConjugate_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeDomainArray.h"
#include"APKCountCells.h"

#include"APKUpdateVarConjugate.h"
#include"APKUpdateVarConjugate2d.h"
#include"APKRefreshSurfaceCells.h"
#include"APKExtractVariables.h"

#include "APKFillNewsolToPrevsol.h"
double CalcVariableError(int variableTag,Base *base);

int CalcConjugate(Base *base,int variableTag);

#endif

