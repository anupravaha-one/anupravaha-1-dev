#ifndef __APKboundary_d_t_sInsulation_h
#define __APKboundary_d_t_sInsulation_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"
#include "APKCalcDiffusionCoefCommon.h"
#include "APKCalGradRad.h"
#include "APKMakeSlabPtrArray.h"
int boundary_d_t_sInsulation(Base *base, int direction, Zone *z, ZoneBoundary *zBc, Dataset *dset, double *dummy);

#endif
