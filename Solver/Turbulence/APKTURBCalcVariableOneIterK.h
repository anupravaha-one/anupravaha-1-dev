#ifndef __APKTURBCalcVariableOneIterK_h
#define __APKTURBCalcVariableOneIterK_h

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
#include "APKTURBOmega.h"
#include "APKTURBOmegaSST.h"
#include "APKUpdateVar.h"
#include "APKUpdateVarCD2d.h"

double CalcVariableOneIterK(Base *base,int variableTag);

#endif



