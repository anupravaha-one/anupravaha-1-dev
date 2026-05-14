#ifndef __APKCalcVariablePotential_h
#define __APKCalcVariablePotential_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"
#include "APKMakeDomainArray.h"
#include "APKCountCells.h"

#include "APKRefreshSurfaceCells.h"
//#include "APKUpdateVarPotential.h"
#include "APKExtractVariables.h"
#include "APKUpdateVarCD2d.h"
#include "APKUpdateVarCD.h"

#include"APK_File_Open.h"

int CalcVariablePotential(Base *base,int variableTag);

#endif
