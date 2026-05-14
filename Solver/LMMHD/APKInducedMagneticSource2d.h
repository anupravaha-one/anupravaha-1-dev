#ifndef __APKInducedMagneticSource2d_h
#define __APKInducedMagneticSource2d_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeDomainArray.h"
#include"APKMakeSlabPtrArray.h"


double InducedMagneticSource2d(Zone *zone,int NI,int NJ,int NK,int index[3],double grad[2]);

#endif

