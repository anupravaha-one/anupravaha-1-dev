#ifndef __APKCalcExtraTermsOfVarDensity2d_h
#define __APKCalcExtraTermsOfVarDensity2d_h

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

int udotdelRho2d(Zone *zone,double *sumU_pdotdelRho);
double CalcExtraTermsOfVarDensity2d(Base *base,Zone *zone);

#endif

