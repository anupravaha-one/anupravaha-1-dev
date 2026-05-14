#ifndef __ElectricFieldXYZ_h
#define __ElectricFieldXYZ_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "../Common/APKData_Structures.h"
#include "../Common/APKMakeSlabPtrArray.h"

#include "APKMakeDomainArray.h"
#include "APKCalcPressureGrad.h"
#include "APKRefreshSurfaceCells.h"


int ElectricFieldXYZ(Base *base);

#endif

