#ifndef __APKMakeSlabPtrArray_h
#define __APKMakeSlabPtrArray_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKFindSlabPtrIface.h"
#include"APKFindSlabPtrBC.h"
#include"APKFindInterfaceBCPtrMod.h"

int MakeSlabPtrArray(Zone *,double *ptrArray[6],int );
int MakeSlabPtrArrayDeltaInverse(Zone * ,double *ptrArray[6]);


#endif

