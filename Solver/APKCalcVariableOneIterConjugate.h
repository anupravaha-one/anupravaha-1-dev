#ifndef __APKCalcVariableOneIterConjugate_h
#define __APKCalcVariableOneIterConjugate_h

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
#include"APKUpdateVarCD.h"
#include"APKUpdateVarCD2d.h"
#include"APKUpdateVarQVK.h"
#include"APKUpdateVarQVK2d.h"
#include"APKUpdateVarConjugate.h"
#include"APKUpdateVarConjugate2d.h"

double CalcVariableOneIterConjugate(Base *base,int variableTag);

#endif

