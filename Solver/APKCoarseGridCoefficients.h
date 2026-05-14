#ifndef __APKCoarseGridCoefficients_h
#define __APKCoarseGridCoefficients_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"

#include"APKMakeDiffusionCoeffDomain.h"

double CoarseGridCoefficients(Base *base,Zone *zone,int g);

#endif

