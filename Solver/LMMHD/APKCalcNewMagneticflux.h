#ifndef __APKCalcNewMagneticflux_h
#define __APKCalcNewMagneticflux_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include "APKMakeDomainArray.h"
#include "APKMakeSlabPtrArray.h"

// #include "APKDiffusion3dWithConstDiffCoeff.h"
// #include "APKDiffusion2dWithConstDiffCoeff.h"
#include "APKDiffusion2dWithConstDiffCoeff_opti.h"
#include "APKDiffusion3dWithConstDiffCoeff_opti.h"
#include "APKMakeDiffusionCoeffDomain.h"

int CalcNewMagneticflux(Base *base);
int CalcNewMagneticflux2d(Base *base);

#endif



