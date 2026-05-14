#ifndef __APKShearStress_h
#define __APKShearStress_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKFindInterfaceBCPtrMod.h"
#include"APKMakeDomainArray.h"
#include"APKMakeSlabPtrArray.h"

#include"APKDiffusion.h"
#include"APKCalcDiffusionCoef.h"

int ShearStress(Base *base,int vartag);

#endif

