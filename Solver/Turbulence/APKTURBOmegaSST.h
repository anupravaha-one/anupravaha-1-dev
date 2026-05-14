#ifndef __APKTURBOmegaSST_h
#define __APKTURBOmegaSST_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"

#include"APKDiffusion2dWithVaryingDiffCoeff.h"
#include"APKDiffusion3dWithVaryingDiffCoeff.h"
#include"APKCalcvarCoeffAndPrevCoeff.h"
#include"APKCalcSource.h"
#include"APKconvfluxQ2D.h"
#include"APKconvfluxQ3D.h"
#include"APKconvfluxC3D.h"
#include"APKconvfluxC2D.h"
#include"APKMakeDiffusionCoeffDomain.h"
#include"APKMakeDiffusionCoeffSlabPtrArray.h"
#include "APKTURBCalcGradSquare.h"
#include"APKTURBChangeIndices.h"

double OmegaSST(Base *base,Zone *zone,double *variable,double *prevTstepVar,int vartag,double *prevIterVar);

#endif

