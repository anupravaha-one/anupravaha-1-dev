#ifndef __APKUpdatePressureVarDensity_h
#define __APKUpdatePressureVarDensity_h

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
#include"APKMatrixMultCorner.h"
#include"APKMatrixMult.h"

double UpdatePressureVarDensity(Base *base,Zone *zone);

#endif

