#ifndef __APKUpdatePotential2d_h
#define __APKUpdatePotential2d_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"

#include"APKDiffusion2dWithConstDiffCoeff.h"
#include"APKTDMA.h"
#include"APKMakeDiffusionCoeffDomain.h"
#include"APKCalcSource.h"

double UpdatePotential2d(Base *base,Zone *zone);

#endif

