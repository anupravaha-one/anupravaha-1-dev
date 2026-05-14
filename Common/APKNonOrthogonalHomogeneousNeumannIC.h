#ifndef __APKNonOrthogonalHomogeneousNeumannIC_h
#define __APKNonOrthogonalHomogeneousNeumannIC_h

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

double NonOrthogonalHomogeneousNeumannIC(Base *base, Zone *z,Interface *facePtr,int vartag);
#endif
