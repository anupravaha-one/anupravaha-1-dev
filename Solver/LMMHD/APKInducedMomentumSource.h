#ifndef __APKInducedMomentumSource_h
#define __APKInducedMomentumSource_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeDomainArray.h"
#include"APKMakeSlabPtrArray.h"

double InducedMomentumSource(Zone *zone,int NI,int NJ,int NK,int index[3],double grad[3]);


#endif

