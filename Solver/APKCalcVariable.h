#ifndef __APKCalcVariable_h
#define __APKCalcVariable_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"
#include "APKMakeDomainArray.h"
#include "APKCountCells.h"

#include "APKfillnewtoprev.h"
#include "APKRefreshSurfaceCells.h"
#include "APKUpdateVar.h"
#include "APKUpdateVarCD.h"
#include "APKUpdateVarCD2d.h"
#include "APKUpdateVarQVK.h"
#include "APKUpdateVarQVK2d.h"
#include "APKExtractVariables.h"
#include "APKupdate_dele_iteration.h"
#include "APKfind_convcri.h"
#include "APKCalcSource_opti.h"
#include "APKLimitVariable.h"


double CalcVariable(Base *base,int variableTag);

#endif

