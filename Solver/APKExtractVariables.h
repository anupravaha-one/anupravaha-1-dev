#ifndef __APKExtractVariables_h
#define __APKExtractVariables_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKupdate_boundary.h"

#include "APKData_Structures.h"
#include "APKFillerConjugateInterface.h"
#include "APKFillerConjugateInterfaceTurb.h"
#include "APKFillerIfaceConjugatePressure.h"
#include "APKFillerIfaceConjugateTemperatureTurb.h"
#include "APKFillerIfaceConjugateVariable.h"
#include "APKFillerIfaceConjugatePotentialMHD.h"
#include "APKFillerInterfaceVariables.h"
#include "APKFillerInterfaceRadiationVariables.h"
#include "APKfind_convcri.h"
#include "APKFindInterfaceBCPtrMod.h"
#include "APKFindZonePtr.h"
#include "APKHomogeneousDirichletIC.h"
// #include "APKHomogeneousNeumannIC.h"
#include "APKNonOrthogonalHomogeneousNeumannIC.h"
#include "APKNonOrthogonalNeumannIC.h"
#include"APKvelocitymassBc.h"

int ExtractVariables(int variableTag,Base *base,Zone *zone);

#endif

