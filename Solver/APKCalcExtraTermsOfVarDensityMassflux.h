#ifndef __APKCalcExtraTermsOfVarDensityMassflux_h
#define __APKCalcExtraTermsOfVarDensityMassflux_h

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
#include"APKCalcExtraTermsOfVarDensityMassflux2d.h"

double CalcExtraTermsOfVarDensityMassflux(Base *base,Zone *zone);

#endif

