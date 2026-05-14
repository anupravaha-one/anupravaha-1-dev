#ifndef __APKFillerIfaceConjugateVariable_h
#define __APKFillerIfaceConjugateVariable_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKCalcDiffusionCoefCommon.h"
#include "APKData_Structures.h"
#include "APKDiffusionCommon.h"
#include "APKfind_convcri.h"
#include "APKFindConjuIfaceNeighbor.h"
#include "APKMakeDomainArray.h"
#include "APKMakeSlabPtrArray.h"
#include "APKTransformCellIndicesCommon.h"

int FillerIfaceConjugateVariable(Base *,int ,Zone *,Zone *,Interface *,int );

#endif

