#ifndef __APKHomogeneousNeumannIC_h
#define __APKHomogeneousNeumannIC_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"
#include "APKFindSlabPtrIface.h"
#include "APKMakeDomainArray.h"


int HomogeneousNeumannIC(Zone *zone,Interface *facePtr,int variableTag);

#endif
