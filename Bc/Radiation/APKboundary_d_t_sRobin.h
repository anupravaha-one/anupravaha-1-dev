#ifndef __APKboundary_d_t_sRobin_h
#define __APKboundary_d_t_sRobin_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"
#include "APKMakeSlabPtrArray.h"
#include "APKCalcDiffusionCoefCommon.h"
#include "APKCalGradRad.h"

int boundary_d_t_sRobin(Base *base , int direction, Zone *z,ZoneBoundary *zBc, Dataset *dset, double *dummy);

#endif
