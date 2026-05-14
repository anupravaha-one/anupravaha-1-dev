#ifndef __APKExtractVariableCoefficients_h
#define __APKExtractVariableCoefficients_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKupdateBCslabsVarCoeff.h"

#include "APKData_Structures.h"
#include "APKFindInterfaceBCPtrMod.h"
#include "APKFillerInterfaceVariables.h"

#include "APKFindZonePtr.h"

int  ExtractVariableCoefficients(int variableTag,Base *base,Zone *zone);

#endif



