#ifndef __APKUpdatePressureArtificial_h
#define __APKUpdatePressureArtificial_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"

#include"APKMakeDiffusionCoeffDomain.h"
#include"APKDiffusion3dWithConstDiffCoeff.h"
#include"APKMatrixMultCorner.h"
#include"APKMatrixMult.h"

double UpdatePressureArtificial(Base *base,Zone *zone);
int CountCells1(Zone *zone);

#endif

