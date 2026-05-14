#ifndef __APKTURBCalcVariableOneIterOmega_h
#define __APKTURBCalcVariableOneIterOmega_h

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
#include "APKTURBOmega.h"
#include "APKUpdateVar.h" 
#include "APKUpdateVarCD2d.h" 

double CalcVariableOneIterOmega(Base *base,int variableTag);

#endif
