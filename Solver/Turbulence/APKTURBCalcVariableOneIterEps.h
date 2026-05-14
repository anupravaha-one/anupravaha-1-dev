#ifndef __APKTURBCalcVariableOneIterEps_h
#define __APKTURBCalcVariableOneIterEps_h

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
#include "APKTURBkturb.h"
#include "APKUpdateVar.h"

double CalcVariableOneIterEps(Base *base,int variableTag);

#endif

