#ifndef __APKCNterm_h
#define __APKCNterm_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"
#include"APKMakeDomainArray.h"

#include"APKfillnewtoprev.h"
#include"APKconvfluxQ2D.h"
#include"APKconvfluxQ3D.h"
#include"APKconvfluxC3D.h"
#include"APKconvfluxC2D.h"
#include"APKDiffusion2dWithConstDiffCoeff.h"
#include"APKDiffusion2dWithVaryingDiffCoeff.h"
#include"APKDiffusion3dWithConstDiffCoeff.h"
#include"APKDiffusion3dWithVaryingDiffCoeff.h"
#include"APKMakeDiffusionCoeffSlabPtrArray.h"
#include"APKMakeDiffusionCoeffFlag.h"
#include"APKMakeDiffusionCoeffDomain.h"

int UpdateCNterm(Base *base,Zone *zone,int vartag);
int CalcPrevTStepFluxes(Base *base);
double *MakeDomainArrayCNTerm(Zone *zone,int vartag);

#endif

