#ifndef __APKCalcVariableCharge_h
#define __APKCalcVariableCharge_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"

#include "APKCountCells.h"

#include "APKRefreshSurfaceCells.h"

#include "APKRefreshWireCharge.h"

#include "APKUpdateVarCD2d.h"
#include "APKUpdateVarCD.h"
#include "APKUpdateVarQVK2d.h"
#include "APKUpdateVarQVK.h"

#include "APKExtractVariables.h"

#include"APK_File_Open.h"

int CalcVariableCharge(Base *base,int variableTag);

#endif

