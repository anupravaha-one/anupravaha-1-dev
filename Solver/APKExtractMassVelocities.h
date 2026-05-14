#ifndef __APKExtractMassVelocities_h
#define __APKExtractMassVelocities_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"
#include "APKFillerConjugateInterfaceTurb.h"
#include "APKFillerIfaceConjugateVariable.h"
#include "APKFillerInterfaceVariables.h"
#include "APKFindZonePtr.h"
#include "APKHomogeneousDirichletIC.h"
#include "APKvelocitymassBc.h"

int ExtractMassVelocities(int variableTag,Base *base,Zone *zone);

#endif

