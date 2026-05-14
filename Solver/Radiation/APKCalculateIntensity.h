#ifndef __APKCalculateIntensity_h
#define __APKCalculateIntensity_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"
#include "APKMakeDomainArray.h"
#include "APKSolution3DConstRefIndex.h"
#include "APKSolution2DConstRefIndex.h"

int CalculateIntensity(Base *base, Zone *zone, int variableTag);

#endif
