#ifndef __APKupdate_boundary_h
#define __APKupdate_boundary_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKboundary_d_t_s.h"
#include "APKboundary_d_t_sTurb.h"

#include "APKFind_Direction.h"
#include "APKData_Structures.h"
#include "APKfind_nv.h"
#include "APKMakeDomainArray.h"
#include "APKFindSlabPtrBC.h"

int update_boundary(Base *base,int nv,Zone *z,ZoneBoundary *ZoneBc);
int mulphase(int,int);

#endif

