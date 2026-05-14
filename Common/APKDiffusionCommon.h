#ifndef __APKDiffusionCommon_h
#define __APKDiffusionCommon_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <math.h>

#include "APKData_Structures.h"

double DiffusionCommon(int index[3],double *,double *,int *,double *slabVariable[6],double *slabVolume[6],double *,double DiffFlux[6]);

#endif

