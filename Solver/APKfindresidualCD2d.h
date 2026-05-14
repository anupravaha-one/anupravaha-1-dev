#ifndef __APKfindresidualCD2d_h
#define __APKfindresidualCD2d_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"

#include"APKCalcSource.h"
#include"APKconvfluxQ2D.h"
#include"APKconvfluxQ3D.h"
#include"APKconvfluxC3D.h"
#include"APKconvfluxC2D.h"
#include"APKDiffusion2dWithConstDiffCoeff.h"
#include"APKDiffusion2dWithVaryingDiffCoeff.h"
#include"APKDiffusion3dWithConstDiffCoeff.h"
#include"APKDiffusion3dWithVaryingDiffCoeff.h"
#include"APKExtractMassVelocities.h"
#include"APKMakeDiffusionCoeffDomain.h"
#include"APKMakeDiffusionCoeffSlabPtrArray.h"
#include"APKMakeDiffusionCoeffFlag.h"
#include"APKMatrixMult.h"
#include"APKvelocitymassBc.h"

double findresidualCD2d(Base *base,Zone *zone,int vartag);

#endif

