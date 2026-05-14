#ifndef __APKCalcCurrentDensity_h
#define __APKCalcCurrentDensity_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKFindgradientPhiMHD.h"
#include"../Geometry/APKface_center.h"
#include"../Geometry/APKFindCellCentre.h"
#include"APKMakeSlabPtrArray.h"
#include"APKMakeDomainArray.h"
#include"APKMakeDiffusionCoeffDomain.h"
#include"APKCalcPressureGrad.h"
#include"APKDiffusion2dWithConstDiffCoeff_opti.h"
#include"APKDiffusion3dWithConstDiffCoeff_opti.h"

double CalcCurrentDensity(Base *base);

#endif

