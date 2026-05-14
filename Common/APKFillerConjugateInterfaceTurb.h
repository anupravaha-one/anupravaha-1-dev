#ifndef __APKFillerConjugateInterfaceTurb_h
#define __APKFillerConjugateInterfaceTurb_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"
#include "APKfind_convcri.h"
#include "APKFindSlabPtrIface.h"
#include "APKHomogeneousDirichletIC.h"
#include "APKHomogeneousNeumannIC.h"
#include "APKNonOrthogonalHomogeneousNeumannIC.h"
#include "APKNonOrthogonalNeumannIC.h"
#include "APKWallFunctionDirichletConjugate.h"

int FillerConjugateInterfaceTurb(Base *base,int variableTag,Interface *facePtr,int iface ,Zone *zone);

#endif
