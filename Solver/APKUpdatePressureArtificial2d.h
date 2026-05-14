#ifndef __APKUpdatePressureArtificial2d_h
#define __APKUpdatePressureArtificial2d_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"

#include"APKDiffusion2dWithConstDiffCoeff.h"
//#include"APKDiffusion2d.h" //Not needed
#include"APKTDMA.h"
#include"APKMakeDiffusionCoeffDomain.h"
#include"APKMatrixMultCorner.h"
#include"APKMatrixMult.h"

double UpdatePressureArtificial2d(Base *base,Zone *zone);

#endif

