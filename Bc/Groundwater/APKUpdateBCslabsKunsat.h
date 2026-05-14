#ifndef __APKUpdateBCslabsKunsat_h
#define __APKUpdateBCslabsKunsat_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>



#include "APKData_Structures.h"
#include "APKFindSlabPtrBC.h"
#include "APKfind_nv.h"

int UpdateBCslabsKunsat(Base *base,Zone *z,ZoneBoundary *ZoneBc,int nv,int vartag,int direct_index);
double	FindKunsat(int indexH,double *Head,Zone *z,ZoneBoundary *ZoneBc);

#endif

