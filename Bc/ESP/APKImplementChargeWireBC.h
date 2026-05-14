#ifndef __APKImplementChargeWireBC_h
#define __APKImplementChargeWireBC_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "../Common/APKData_Structures.h"
#include "../Common/APKMakeSlabPtrArray.h"
#include "../Common/APKMakeDomainArray.h"
#include "APKDiffusion2d.h"

int ImplementChargeWireBC(Base* base ,int direction,double *flowsol,int vartag,Zone *z,int size,Dataset *dset,double *dummy);

#endif

