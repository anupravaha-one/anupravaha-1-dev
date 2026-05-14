#ifndef __APKFillerIfaceConjugateTemperatureTurb_h
#define __APKFillerIfaceConjugateTemperatureTurb_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKCalcDiffusionCoefCommon.h"
#include "APKData_Structures.h"
#include "APKDiffusionCommon.h"
#include "APKDuplicateSlabTemperatures.h"
#include "APKfind_convcri.h"
#include "APKFindConjuIfaceNeighbor.h"
#include "APKFindInterfaceBCPtrMod.h"
#include "APKMakeDomainArray.h"
#include "APKMakeSlabPtrArray.h"
#include "APKModifiedDiffCoeff.h"
#include "APKTransformCellIndicesCommon.h"
#include "APKWallFunctionDirichletConjugate.h"

int FillerIfaceConjugateTemperatureTurb(Base *,int ,Zone *,Zone *,Interface *);

#endif

