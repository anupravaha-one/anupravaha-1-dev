#ifndef __APKupdate_boundaryTurbDMPFV_h
#define __APKupdate_boundaryTurbDMPFV_h


#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKfind_nv.h"
#include "../Common/APKData_Structures.h"

int update_boundaryTurbDMPFV(Base *base,Zone *CgnsZone,ZoneBoundary *bcptr,int vartag);

#endif
