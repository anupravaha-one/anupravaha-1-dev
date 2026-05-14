#ifndef __APKCalcExtraTermsOfVarDensityMassflux2d_h
#define __APKCalcExtraTermsOfVarDensityMassflux2d_h

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

double CalcExtraTermsOfVarDensityMassflux2d(Base *base,Zone *zone);
int CalcMassFlux(double *newMassflux, Zone *zone);

#endif

