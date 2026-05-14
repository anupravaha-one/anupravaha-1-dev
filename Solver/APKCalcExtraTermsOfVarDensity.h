#ifndef __APKCalcExtraTermsOfVarDensity_h
#define __APKCalcExtraTermsOfVarDensity_h

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

double CalcExtraTermsOfVarDensity(Base *base,Zone *zone);
int udotdelRho(Zone *zone,double *sumU_pdotdelRho);

#endif

