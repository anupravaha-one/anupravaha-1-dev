#ifndef __APKSolution3DConstRefIndex_h
#define __APKSolution3DConstRefIndex_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"
#include "APKMakeSlabPtrArray.h"
#include "APKChoiceScatteringfunction.h"
#include "APKRadiationSolutionInt3D.h"
#include "APKIntExtremeFace3D.h"
#include "APKIntExtremeCell3D.h"

int Solution3DConstRefIndex(Base *base, Zone *zone, double *variableDomain, double *prevVariableDomain, double alpha);

#endif
