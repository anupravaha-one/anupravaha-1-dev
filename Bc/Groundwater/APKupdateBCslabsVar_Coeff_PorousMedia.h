#ifndef __APKupdateBCslabsVar_Coeff_PorousMedia_h
#define __APKupdateBCslabsVar_Coeff_PorousMedia_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"
#include "APKFindSlabPtrBC.h"

#include "APKGetSurfacePoints.h"
#include "APKFindFaceCenterBC.h"
#include "APKfind_nv.h"


int updateBCslabsVar_Coeff_PorousMedia(Base *base,int variableTag,Zone *zone,ZoneBoundary *bcPtr,int dir);
int UpdateBcDarcyVelocity(Base *base,int variableTag,Zone *zone,ZoneBoundary *bcPtr,int dir);
int updateBcDispersionCoeff(Base *base,int variableTag,Zone *zone,ZoneBoundary *bcPtr,int dir);

#endif

