#ifndef __APKRadiationSolverWithFlow3D_h
#define __APKRadiationSolverWithFlow3D_h

#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>
#include<string.h>

#include "APKData_Structures.h"
#include "APKFillNewRadsolToPrevRadsol.h"
#include "APKCalculateIntensity.h"
#include "APKFillOutgoingZoneBCIntensity.h"
#include "APKExtractVariables.h"

int RadiationSolverWithFlow3D(Base *CgnsBase);

#endif
