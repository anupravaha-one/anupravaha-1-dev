#ifndef __APKUpdatePotentialMHD_h
#define __APKUpdatePotentialMHD_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"

#include"APKDiffusion3dWithConstDiffCoeff.h"
#include"APKMakeDiffusionCoeffDomain.h"

double UpdatePotentialMHD(Base *base,Zone *zone);

#endif

