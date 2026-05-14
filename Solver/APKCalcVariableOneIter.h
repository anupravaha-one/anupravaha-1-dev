#ifndef __APKCalcVariableOneIter_h
#define __APKCalcVariableOneIter_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeDomainArray.h"

#include"APKExtractVariables.h"
#include"APKfillnewtoprev.h"
#include"APKUpdateVar.h"
#include"APKUpdateVarQVK.h"
#include"APKUpdateVarQVK2d.h"
#include"APKUpdateVarCD.h"
#include"APKUpdateVarCD2d.h"
#include"APKRefreshSurfaceCells.h"
#include"APKUpdateScaler.h"

double CalcVariableOneIter(Base *base,int variableTag);
double CalcVariableOneIter_porousMedia(Base *base,int variableTag);

#endif

