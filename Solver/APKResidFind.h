#ifndef __APKResidFind_h
#define __APKResidFind_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APKMakeSlabPtrArray.h"

#include"APKMatrixMult.h"
#include"APKMatrixMultCorner.h"

double ResidFind(Base *base,Zone *zone,int g);

#endif

