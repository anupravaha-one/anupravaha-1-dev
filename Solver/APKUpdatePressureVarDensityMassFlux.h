#ifndef __APKUpdatePressureVarDensityMassFlux_h
#define __APKUpdatePressureVarDensityMassFlux_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"
#include"APKMakeDomainArray.h"

#include"APKuserDefinedVarCoefficients.h"

#include"APKDiffusion3dWithConstDiffCoeff.h"
#include"APKMakeDiffusionCoeffDomain.h"
#include"APKUpdatePressure2dVarDensityMassFlux.h"
#include"APKMatrixMultCorner.h"
#include"APKMatrixMult.h"

double UpdatePressureVarDensityMassFlux(Base *base,Zone *zone);


#endif

