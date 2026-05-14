#ifndef __APKMakeDiffusionCoeffSlabPtrArray_h
#define __APKMakeDiffusionCoeffSlabPtrArray_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKFindInterfaceBCPtrMod.h"

#include"APKFindDiffusionCoeffSlabPtrIface.h"
#include"APKFindDiffusionCoeffSlabPtrBC.h"

int MakeDiffusionCoeffSlabPtrArray(Zone *zone,double *ptrArray[6],int vartag);

#endif

