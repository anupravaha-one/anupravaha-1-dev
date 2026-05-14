#ifndef __APKvelocitymassBc_h
#define __APKvelocitymassBc_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKFind_Direction.h"

#include"APKboundary_d_t_s.h"
#include"APKboundary_d_t_sTurb.h"
#include"APKfind_nv.h"
#include"APKupdate_boundary.h"

int velocitymassBc(Base *base,int nv,Zone *z,ZoneBoundary *ZoneBc);

#endif

