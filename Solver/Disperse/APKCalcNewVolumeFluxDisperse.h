#ifndef __APKCalcNewVolumeFluxDisperse_h
#define __APKCalcNewVolumeFluxDisperse_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"

#include "APKMakeDomainArray.h"
#include "APKMakeSlabPtrArray.h"

#include "APKDiffusion3dWithConstDiffCoeff.h"
#include "APKDiffusion2dWithConstDiffCoeff.h"

#include "APKDiffusion2dWithConstDiffCoeff_opti.h"
#include "APKDiffusion3dWithConstDiffCoeff_opti.h"

#include "APKDiffusion2dWithConstDiffCoeff_optiMPFV.h"
#include "APKDiffusion3dWithConstDiffCoeff_optiMPFV.h"

#include "APKMakeDiffusionCoeffDomain.h"

int CalcNewVolumeFluxDisperse(Base *);
int CalcNewVolumeFluxDisperse2d(Base *);
#endif
