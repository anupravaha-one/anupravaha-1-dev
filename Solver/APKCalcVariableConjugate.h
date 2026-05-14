#ifndef __APKCalcVariableConjugate_h
#define __APKCalcVariableConjugate_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeDomainArray.h"
#include"APKCountCells.h"

#include"APKfillnewtoprev.h"
#include "APKfind_convcri.h"
#include"APKUpdateVarConjugate.h"
#include"APKUpdateVarConjugate2d.h"
#include"APKRefreshSurfaceCells.h"
#include"APKUpdateVar.h"
#include"APKUpdateVarCD.h"
#include"APKUpdateVarCD2d.h"
#include"APKUpdateVarQVK.h"
#include"APKUpdateVarQVK2d.h"
#include"APKExtractVariables.h"

int CalcVariableConjugate(Base *base,int variableTag);

#endif

