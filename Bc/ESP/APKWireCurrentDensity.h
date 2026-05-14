#ifndef __WireCurrentDensity_h
#define __WireCurrentDensity_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "../Common/APKData_Structures.h"
#include "../Common/APKMakeSlabPtrArray.h"
#include "../Common/APKCalcDiffusionCoefCommon.h"
#include "../Common/APKFindSlabPtrBC.h"

#include "APKFindInterfaceBCPtrMod.h"
#include "APKMakeDomainArray.h"
#include "APKDiffusion.h"
#include "APKDiffusion2d.h"

#include "../Geometry/APKFindCellCentre.h"

int WireCurrentDensity(Base *base,int vartag);

#endif

