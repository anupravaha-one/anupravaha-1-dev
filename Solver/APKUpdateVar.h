#ifndef __APKUpdateVar_h
#define __APKUpdateVar_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"

#include"APKMakeSlabPtrArray.h"

#include"APKDiffusion2dWithConstDiffCoeff.h"
#include"APKDiffusion2dWithVaryingDiffCoeff.h"
#include"APKDiffusion3dWithConstDiffCoeff.h"
#include"APKDiffusion3dWithVaryingDiffCoeff.h"
#include"APKCalcvarCoeffAndPrevCoeff.h"
#include"APKCalcSource.h"
#include"APKconvfluxQ2D.h"
#include"APKconvfluxQ3D.h"
#include"APKconvfluxC3D.h"
#include"APKconvfluxC2D.h"
#include"APKMakeDiffusionCoeffSlabPtrArray.h"
#include"APKMakeDiffusionCoeffDomain.h"
#include"APKMakeDiffusionCoeffFlag.h"
#include"APKTURBCalcGradSquare.h"

double UpdateVar(Base *base,Zone *zone,double *variable,double *prevTstepVar,int vartag,double *prevIterVar);

#endif
