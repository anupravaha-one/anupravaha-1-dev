#ifndef __APKTURBgrad2D_h
#define __APKTURBgrad2D_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"
#include "APKMakeSlabPtrArray.h"
#include "APKMakeDomainArray.h"

double grad2D(Base *base,Zone *zone,double *variable,double *slabVariable[6],double *variable1,
		double *slabVariable1[6],int INDEX[3]);


#endif

