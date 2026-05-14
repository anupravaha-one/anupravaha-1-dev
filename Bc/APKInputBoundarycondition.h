#ifndef __APKInputBoundarycondition_h
#define __APKInputBoundarycondition_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"
 
#include "APKcheck_function.h"
#include "APKcheck_bc_type.h"
#include "APKprint_tagname.h"
#include "APKPoint_to_Coord.h"
#include "APKUser_Enter.h"
#include "APKCall_user_func.h"
#include "APKfind_name.h"
#include  "APKupdate_boundaryDMPFV.h"
#include  "APKupdate_boundaryDisperse.h"
#include  "APKupdate_boundary.h"
#include "APKfind_nv.h"

#include "APKboundary_d_t_s.h"
#include "APKFind_Direction.h"
#include "APKMakeDomainArray.h"
#include "APKFindSlabPtrBC.h"


int InputBoundarycondition(Base *base,int NS,Zone *CgnsZone,int flag2d);


#endif

