#ifndef __APKSolution2DConstRefIndex_h
#define __APKSolution2DConstRefIndex_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"
#include "APKMakeSlabPtrArray.h"
#include "APKChoiceScatteringfunction.h"
#include "APKRadiationSolutionInt2D.h"
#include "APKIntExtremeFace2D.h"
#include "APKIntExtremeCell2D.h"

int Solution2DConstRefIndex(Base *base, Zone *zone, double *variableDomain, double *prevVariableDomain, double alpha);

#endif
