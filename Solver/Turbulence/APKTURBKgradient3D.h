#ifndef __APKTURBKgradient3D_h
#define __APKTURBKgradient3D_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"../Common/APKData_Structures.h"
#include"../Common/APKMakeSlabPtrArray.h"
#include"../Common/APKMakeDomainArray.h"

double Kgradient3D(Base *base,Zone *zone,double *variable,double *slabVariable[6],int INDEX[3]);

#endif
