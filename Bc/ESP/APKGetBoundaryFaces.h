#ifndef __APKGetBoundaryFaces_h
#define __APKGetBoundaryFaces_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>


#include "APKData_Structures.h"
#include "APKFindInterfaceBCPtrMod.h"
#include "APKUpdate_Boundary_Charge.h"
#include "APKFindZonePtr.h"

int GetBoundaryFaces(int variableTag,Base *base,Zone *zone);

#endif

