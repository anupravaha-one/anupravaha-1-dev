#ifndef __APKupdateBCDeltaInverse_h
#define __APKupdateBCDeltaInverse_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"
#include "APKTransformCellIndicesCommon.h"


int updateBCDeltaInverse(int c, Zone *zone,ZoneBoundary *bcPtr,int dir);
int FillerInterfaceDelta(int c,Interface *facePtr,Zone *ptrZone,int iface);
int findInterfaceZoneDelta(int c,Interface *facePtr,Zone *zone,int index2[3],int indexofSlabs);

#endif











