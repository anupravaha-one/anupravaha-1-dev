#ifndef __APKupdateBCslabsVarCoeff_h
#define __APKupdateBCslabsVarCoeff_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "../Common/APKData_Structures.h"
#include "../Common/APKFindSlabPtrBC.h"

#include "../Geometry/APKGetSurfacePoints.h"
#include "../Geometry/APKFindFaceCenterBC.h"

#include "../User/APKBC_userDefinedVarCoefficients.h"

int updateBCslabsVarCoeff(int variableTag,Zone *zone,ZoneBoundary *bcPtr,int dir,Base *base);

#endif

