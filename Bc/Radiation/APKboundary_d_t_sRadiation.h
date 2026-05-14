#ifndef __APKboundary_d_t_sRadiation_h
#define __APKboundary_d_t_sRadiation_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"
#include "APKChoiceScatteringfunction.h"
#include "APKMakeDomainArray.h"

int boundary_d_t_sRadiation(Base *base, int direction, Zone *z, ZoneBoundary *zBc, Dataset *dset, double *dummy);

#endif
