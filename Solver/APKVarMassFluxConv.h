#ifndef __APKVarMassFluxConv_h
#define __APKVarMassFluxConv_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKFindInterfaceBCPtrMod.h"
#include"APKMakeDomainArray.h"
#include"APKMakeSlabPtrArray.h"

#include"APKconvfluxC2D.h"
#include"APKDiffusion.h"
#include"APKCalcDiffusionCoef.h"

double VarMassFluxConv(Base *base,int vartag);

#endif

