#ifndef __APKMultigridExtractVariables_h
#define __APKMultigridExtractVariables_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKupdate_boundary.h"

#include "APKData_Structures.h"
#include "APKFindInterfaceBCPtrMod.h"
#include "APKFillerConjugateInterface.h"
#include "APKFillerIfaceConjugatePressure.h"
#include "APKFillerIfaceConjugateVariable.h"
#include "APKFillerInterfaceVariables.h"

#include "APKFindZonePtr.h"

int MultigridExtractVariables(int variableTag,Base *base,Zone *zone,int gdi);

#endif

