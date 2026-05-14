#ifndef __APKUpdatePressureHead2d_h
#define __APKUpdatePressureHead2d_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"

#include"APKDiffusion2dWithConstDiffCoeff.h"
#include"APKDiffusion2dWithVaryingDiffCoeff.h"
#include"APKdatum2d.h"
#include"APKdatumUnsat2d.h"
#include"APKMakeDiffusionCoeffSlabPtrArray.h"
#include"APKMakeDiffusionCoeffDomain.h"


double UpdatePressureHead2d(Base *base,Zone *zone);

#endif
