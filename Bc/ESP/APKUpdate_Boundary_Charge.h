#ifndef __APKUpdate_Boundary_Charge_h
#define __APKUpdate_Boundary_Charge_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>


#include "APKImplementChargeWireBC.h"
#include "APKWriteWireCharge.h"

#include "../Common/APKFind_Direction.h"
#include "../Common/APKData_Structures.h"
#include "../Common/APKfind_nv.h"

int Update_Boundary_Charge(Base *base,int nv,Zone *z,ZoneBoundary *ZoneBc);

#endif

