#ifndef __APKFillerInterfaceRadiationVariables_h
#define __APKFillerInterfaceRadiationVariables_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"
#include "APKTransformCellIndicesCommon.h"
#include "APKFindSlabPtrIface.h"
#include "APKMakeDomainArray.h"

int FillerInterfaceRadiationVariables(Base *base, int variableTag,Interface *facePtr,Zone *ptrZone,int iface);
#endif
