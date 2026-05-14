#ifndef __APKUpdatePressureHead3D_h
#define __APKUpdatePressureHead3D_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"




#include"APKdatum2d.h"
#include"APKdatumUnsat2d.h"
#include"APKMakeDiffusionCoeffSlabPtrArray.h"
#include"APKMakeDiffusionCoeffDomain.h"
#include"APKDiffusion3dWithConstDiffCoeff.h"   
#include"APKDiffusion3dWithVaryingDiffCoeff.h"
  

double UpdatePressureHead3D(Base *base,Zone *zone);

#endif
