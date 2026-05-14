#ifndef __APKUpdatePressure2dVarDensity_h
#define __APKUpdatePressure2dVarDensity_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"
#include"APKMakeDomainArray.h"

#include"APKuserDefinedVarCoefficients.h"

#include"APKDiffusion2dWithConstDiffCoeff.h"
#include"APKMakeDiffusionCoeffDomain.h"
#include"APKMatrixMultCorner.h"
#include"APKMatrixMult.h"

double UpdatePressure2dVarDensity(Base *base,Zone *zone);

#endif

