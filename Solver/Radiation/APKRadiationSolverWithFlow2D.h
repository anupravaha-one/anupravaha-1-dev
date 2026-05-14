#ifndef __APKRadiationSolverWithFlow2D_h
#define __APKRadiationSolverWithFlow2D_h

#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>
#include<string.h>

#include "APKData_Structures.h"
#include "APKFillNewRadsolToPrevRadsol.h"
#include "APKCalculateIntensity.h"
#include "APKFillOutgoingZoneBCIntensity.h"
#include "APKExtractVariables.h"

int RadiationSolverWithFlow2D(Base *CgnsBase);

#endif
