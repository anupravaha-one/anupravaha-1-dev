#ifndef __APKUpdatePressure2d_h
#define __APKUpdatePressure2d_h

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
#include"APKMatrixMult.h"
#include"APKMatrixMultCorner.h"
#include"APKFindInterfaceBCPtrMod.h"

double UpdatePressure2d(Base *base,Zone *zone);

#endif

