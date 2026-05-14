#ifndef __APKTURBCalcVariableOneIterOmegaSST_h
#define __APKTURBCalcVariableOneIterOmegaSST_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"
#include "APKMakeDomainArray.h"

#include "APKExtractVariables.h"
#include "APKfillnewtoprev.h"
#include "APKRefreshSurfaceCells.h"
#include "APKTURBOmegaSST.h"
#include "APKUpdateVar.h"

double CalcVariableOneIterOmegaSST(Base *base,int variableTag);

#endif
