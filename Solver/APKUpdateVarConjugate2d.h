#ifndef __APKUpdateVarConjugate2d_h
#define __APKUpdateVarConjugate2d_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"

#include"APKCalcvarCoeffAndPrevCoeffConjugate.h"
#include"APKCalcSource.h"
#include"APKDiffusion2dWithConstDiffCoeff.h"
#include"APKDiffusion2dWithVaryingDiffCoeff.h"
#include"APKDiffusion3dWithConstDiffCoeff.h"
#include"APKDiffusion3dWithVaryingDiffCoeff.h"
#include"APKMakeDiffusionCoeffSlabPtrArray.h"
#include"APKMakeDiffusionCoeffDomain.h"
#include"APKMakeDiffusionCoeffFlag.h"
#include"APKMatrixMultCorner.h"
#include"APKMatrixMult.h"

double UpdateVarConjugate2d(Base *base,Zone *zone,int vartag);

#endif

