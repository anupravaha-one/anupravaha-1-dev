#ifndef __APKInterZoneVariableTransfer_h
#define __APKInterZoneVariableTransfer_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"
#include "APKFindInterfaceBCPtrMod.h"
#include "APKFillerConjugateInterface.h"
#include "APKFillerIfaceConjugatePressure.h"
#include "APKFillerIfaceConjugateVariable.h"
#include "APKFillerInterfaceVariables.h"
#include "APKMFillerInterfaceVariables.h"

#include "APKFindZonePtr.h"

int InterZoneVariableTransfer(int variableTag,Base *base,Zone *zone,int gdi);

#endif

