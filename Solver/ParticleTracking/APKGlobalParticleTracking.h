#ifndef __APKGlobalParticleTracking_h
#define __APKGlobalParticleTracking_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
//#include"APKMakeDomainArray.h"
#include"APKFindZonePtr.h"
#include"APKFindInterfaceBCPtrMod.h"
#include"APKTransformCellIndices1.h"

int GlobalParticleTracking(Base *base,double XP,double YP,double ZP);

#endif
