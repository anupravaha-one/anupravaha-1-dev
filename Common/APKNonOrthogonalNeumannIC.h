#ifndef __APKNonOrthogonalNeumannIC_h
#define __APKNonOrthogonalNeumannIC_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKCalcDiffusionCoefCommon.h"
#include "APKData_Structures.h"
#include "APKDiffusionCommon.h"
#include "APKfind_convcri.h"
#include "APKFindSlabPtrIface.h"
#include "APKMakeDomainArray.h"
#include "APKMakeSlabPtrArray.h"

double NonOrthogonalNeumannIC(Base *base, Zone *z,Interface *facePtr,int vartag, double inwardFlux);
#endif
