#ifndef __APKRadiationSolverFiniteVolume2D_h
#define __APKRadiationSolverFiniteVolume2D_h

#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>
#include<string.h>

#include "APKData_Structures.h"
#include "APKFillNewRadsolToPrevRadsol.h"
#include "APKCalculateIntensity.h"
#include "APKFillOutgoingZoneBCIntensity.h"
#include "APKExtractVariables.h"

int RadiationSolverFiniteVolume2D(Root CgnsRoot);

#endif
